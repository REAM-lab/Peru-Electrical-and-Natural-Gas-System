println(string(""))
println(string("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"))
println(string("Energy Integration SP y RGN Acoplado"))
println(string("Modelo eMISOCP"))
println(string("by Pablo Salas"))

#Packages
using JuMP,LinearAlgebra,DataFrames,CSV,Dates,Gurobi;

ss(dataframe, header, name, column)= dataframe[dataframe[:,header].== name, column][1]

#Load directory path
paths=CSV.read("directory_paths.csv",DataFrame);
folder_name=pwd();
#folder_name=paths[1,2];
input_folder=paths[2,2];
output_folder=paths[3,2];

#presion unit, formulación con kPa cambian Cmn y Kmn # Cambiar para caso SEIN_Estudio
cons_rgn_unit_input="kpa"

#condicones iniciales ini or fin
con_ini="fin"
lpini=1700

#Input data
path_inputdata=string(folder_name, "/",input_folder,"/");

#Eficiencia
cons_efi="Eficiencia"

#Hard Code
folder_end="Acoplado_eMISOCP"
#Output data folder
path_output=string(folder_name, "/",output_folder, "/",folder_end,"/");

#Load Configuration.csv
config=CSV.read(string(path_inputdata,"configuration.csv"), DataFrame);
case_name=config[1,2];
start_case=config[2,2];
if con_ini=="ini"    
    final_case=config[2,2];
else
    final_case=config[3,2];
end
date_format=config[4,2];
resolution=parse(Int,config[5,2]);
tolerance_solver=parse(Float64,config[8,2]);
cores_solver=parse(Float64, config[9,2]);
MVAbase=parse(Float64, config[10,2]);
digits_round=parse(Int, config[11,2]);
price_shed=parse(Int, config[12,2]);
g_price_shed=parse(Int,config[13,2])*1000;
thr=parse(Int,config[14,2]);

#Set main parameters
start_case=Date(start_case,date_format);
final_case=Date(final_case,date_format);
ndays=Dates.value(final_case-start_case)+1;
deltaT=resolution/60;
rel_h=resolution*60/(10^6);
#= stageG=24/deltaT =#
stageG=24;
nper=Int(24/deltaT);

#Upload plant.csv
plant=CSV.read(string(path_inputdata,"ps_ugs.csv"), DataFrame);
ngen=size(plant,1);
insertcols!(plant, 1, :order => 1:ngen);

#Ciclo combinado
ciclo=plant[:, [:CCiclo]];
ciclo=unique(ciclo);
ciclo=ciclo[ciclo.CCiclo .!= "0", :];
ncc=size(ciclo,1);
insertcols!(ciclo, 1, :order => 1:ncc);

#Upload plant hidro
planth=CSV.read(string(path_inputdata,"ps_ugh.csv"), DataFrame);
ngh=size(planth,1);
insertcols!(planth, 1, :order => 1:ngh);
hipas=planth[ planth.emb .== 0, :];
hiemb=planth[ planth.emb .== 1 , :];

#Create a time chart.csv to know the time horizon and duration of each period
time_chart=Array{Any}(nothing,ndays*nper+1,8);
time_chart[1,1:8]=["nday" "nper" "year" "month" "day" "hour" "min" "date"];

    for d in 1:ndays, t in 1:nper    
            time_chart[t+1+(d-1)*nper,1]=d;
            time_chart[t+1+(d-1)*nper,2]=(d-1)*24+t;
            time_chart[t+1+(d-1)*nper,3]=Dates.year(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,4]=Dates.month(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,5]=Dates.day(start_case+Dates.Day(d-1));
            time_chart[t+1+(d-1)*nper,6]=Dates.hour(DateTime(start_case)+Dates.Minute(resolution*(t-1)));
            time_chart[t+1+(d-1)*nper,7]=Dates.minute(DateTime(start_case)+Dates.Minute(resolution*(t-1)));        
            time_chart[t+1+(d-1)*nper,8]=DateTime(start_case+Dates.Day(d-1))+Dates.Minute(resolution*(t-1));
    end
    time_chart_df = DataFrame(time_chart[2:end,:], :auto);
    rename!(time_chart_df, Symbol.(time_chart[1,:]));
    CSV.write(string(path_output,"time_chart.csv"), time_chart_df);

    time_chart_df[!,:nday] = convert.(Int, time_chart_df[!,:nday]);
    time_chart_df[!,:nper] = convert.(Int, time_chart_df[!,:nper]);
    time_chart_df[!,:year] = convert.(Int, time_chart_df[!,:year]);
    time_chart_df[!,:month] = convert.(Int, time_chart_df[!,:month]);
    time_chart_df[!,:day] = convert.(Int, time_chart_df[!,:day]);
    time_chart_df[!,:hour] = convert.(Int, time_chart_df[!,:hour] );
    time_chart_df[!,:min] = convert.(Int, time_chart_df[!,:min]);

nper=nper*ndays

#Upload buses
bus=CSV.read(string(path_inputdata,"ps_load.csv"), DataFrame);
nbus=size(bus,1);
insertcols!(bus, 1, :order => 1:nbus);
bus
#Upload rer
rer=CSV.read(string(path_inputdata,"ps_rer.csv"), DataFrame)
nrer=size(rer,1);
insertcols!(rer, 1, :order => 1:nrer);
rer
#Load demand
load=CSV.read(string(path_inputdata,"ps_demand.csv"), DataFrame);
load.Date = Date.(load.Date, "dd/mm/yyyy");

#Upload branch
branch=CSV.read(string(path_inputdata,"ps_lines.csv"), DataFrame);
nbranch=size(branch,1);
insertcols!(branch, 1, :order => 1:nbranch);
mapbranch=zeros(nbranch,2);
for br in 1:nbranch
    mapbranch[br,1:2]=[bus[bus.node.==branch[br,:from],:order][1] bus[bus.node.==branch[br,:to],:order][1]];
end
mapbranch=convert.(Int, mapbranch);

#Load G Nodes
gnode=CSV.read(string(path_inputdata,"g_load.csv"), DataFrame);
ngnode=size(gnode,1);
insertcols!(gnode, 1, :order => 1:ngnode);
## Convertir de bar a Kpa
if cons_rgn_unit_input == "bar"
    fac_convert = 1
else
    fac_convert = 100
end

gnode[:,:Pre_max]=gnode[:,:Pre_max]*fac_convert
gnode[:,:Pre_min]=gnode[:,:Pre_min]*fac_convert

#Upload pipes
pipe=CSV.read(string(path_inputdata,"g_pipes.csv"), DataFrame);
npipe=size(pipe,1);
insertcols!(pipe, 1, :order => 1:npipe);
li_pipe=pipe[pipe.Linepack .==1,:]
nli_pipe=pipe[pipe.Linepack .==0,:]

mappipe=zeros(npipe,2);
for pp in 1:npipe
    mappipe[pp,1:2]=[gnode[gnode.Node.==pipe[pp,:First],:order][1] gnode[gnode.Node.==pipe[pp,:Terminal],:order][1]];
end
mappipe=convert.(Int, mappipe);

#Upload wells
well=CSV.read(string(path_inputdata,"g_wells.csv"), DataFrame);
nwell=size(well,1);
insertcols!(well, 1, :order => 1:nwell);
well[:,:C_sup_D]=well[:,:C_sup_D]*1000

#Load demand gas
gload=CSV.read(string(path_inputdata,"g_demand.csv"), DataFrame);
gload.Date = Date.(gload.Date,  "dd/mm/yyyy");

#Profile RER
rer_profile=CSV.read(string(path_inputdata,"ps_profile_rer.csv"), DataFrame);
rer_profile.Date = Date.(rer_profile.Date, "dd/mm/yyyy");

#Cap_Well
capw=CSV.read(string(path_inputdata,"well_cap.csv"), DataFrame);
capw.Date = Date.(capw.Date, "dd/mm/yyyy");

#Upload compressors
comp=CSV.read(string(path_inputdata,"g_compresor.csv"), DataFrame);
ncomp=size(comp,1);
insertcols!(comp, 1, :order => 1:ncomp);
mapcom=zeros(ncomp,2);
for com in 1:ncomp
    mapcom[com,1:2]=[gnode[gnode.Node.==comp[com,:First],:order][1] gnode[gnode.Node.==comp[com,:Terminal],:order][1]];
end
mapcom=convert.(Int, mapcom);

#Upload reg
reg=CSV.read(string(path_inputdata,"g_reg.csv"), DataFrame);
## Convertir de bar a Kpa
#reg[:,:Pin_min]=reg[:,:Pin_min]*fac_convert
#reg[:,:Pin_max]=reg[:,:Pin_max]*fac_convert
#reg[:,:Pout]=reg[:,:Pout]*fac_convert
nreg=size(reg,1);
insertcols!(reg, 1, :order => 1:nreg);
mapreg=zeros(nreg,2);
for re in 1:nreg
    mapreg[re,1:2]=[gnode[gnode.Node.==reg[re,:First],:order][1] gnode[gnode.Node.==reg[re,:Terminal],:order][1]];
end
mapreg=convert.(Int, mapreg);

#Main sets of the optimization problem
NPER=1:nper;
NGEN=1:ngen;
NBUS=1:nbus;
NBR=1:nbranch;
NGN=1:ngnode;
NPP=1:npipe;
NPPL=li_pipe.order;
NPPNL=nli_pipe.order;
NWE=1:nwell;
NCO=1:ncomp;
NRER=1:nrer;
NGH=1:ngh;
NGHP=hipas.order;
NGHE=hiemb.order;
NCC=1:ncc;
NRE=1:nreg;

#Function that finds the active power demand (MW) of the bus n on day d at period t
function pc(n,t)  
    global column=0;
    for i in 1:size(names(load),1)
        if names(load)[i]==string(bus[n, :node])
           global column = i;
           break
        end
   end
   if column !=0 
    d=div(t-1,24)+1
    value=load[ (load.Date .== Date(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])) .& 
    (load.Time .== Time(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])),column][1,1];    
   else
    value=0
   end

    if ismissing(value)
            return 0
    else
            return value
    end
end 

#Function that finds the gas demand of the gnode n on day d at period t
function gc(n,t)  
    global column=0;
    for i in 1:size(names(gload),1)
        if names(gload)[i]==string(gnode[n, :Node])
           global column = i;
           break
        end
   end
   if column !=0 
    d=div(t-1,24)+1
    value=gload[ (gload.Date .== Date(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])) .& 
    (gload.Time .== Time(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])),column][1,1];    
   else
    value=0
   end

    if ismissing(value)
            return 0
    else
            return value
    end
end

#Function that finds the prodution rer on day d at period t
function gr(n,t)  
    global column=0;
    for i in 1:size(names(rer_profile),1)
        if names(rer_profile)[i]==string(rer[n, :unit])
           global column = i;
           break
        end
   end

   if column !=0
        d=div(t-1,24)+1
        value=rer_profile[ (rer_profile.Date .== Date(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])) .& 
        (rer_profile.Time .== Time(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])),column][1,1];    
   else
    value=0
   end 

    if ismissing(value)
            return 0
    else
            return value
    end
end

#Function that finds the prodution well cap on day d at period t
function wc(n,t)  
    global column=0;
    for i in 1:size(names(capw),1)
        if names(capw)[i]==string(well[n, :Well])
           global column = i;
           break
        end
   end

   if column !=0
        d=div(t-1,24)+1
        value=capw[ (capw.Date .== Date(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])) .& 
        (capw.Time .== Time(time_chart_df[ (time_chart_df.nday.==d) .& (time_chart_df.nper.==t), :date][1,1])),column][1,1];    
   else
    value=0
   end 

    if ismissing(value)
            return 0
    else
            return value
    end
end


#build rer matrix
global gerer=zeros(nrer,nper);

loadrer_results=Array{Any}(nothing,nper+2,nrer+1);
loadrer_results[1:2,1]=["order" "rer"];

for t in NPER
 
    loadrer_results[t+2,1]=t;

    for n in NRER
        loadrer_results[1:2,n+1]=[rer[n, :order][1] rer[n, :unit][1]];
     
        gerer[n,t]=gr(n,t);

        loadrer_results[t+2,n+1]=gerer[n,t];              
    
     end
end

#build Capacity well matrix
global well_cap=zeros(nwell,nper);

loadwell_results=Array{Any}(nothing,nper+2,nwell+1);
loadwell_results[1:2,1]=["order" "well"];

for t in NPER
 
    loadwell_results[t+2,1]=t;

    for n in NWE
        loadwell_results[1:2,n+1]=[well[n, :order][1] well[n, :Well][1]];
     
        well_cap[n,t]=wc(n,t);

        loadwell_results[t+2,n+1]=well_cap[n,t];              
    
     end
end

#build power demand matrix
global pcc=zeros(nbus,nper);

loadp_results=Array{Any}(nothing,nper+2,nbus+1);
loadp_results[1:2,1]=["order" "node"];


for t in NPER
 
    loadp_results[t+2,1]=t;

    for n in NBUS
        loadp_results[1:2,n+1]=[bus[n, :order][1] bus[n, :node][1]];
     
        pcc[n,t]=pc(n,t);

        loadp_results[t+2,n+1]=pcc[n,t];              
    
     end
end

#build gas demand matrix

global fl=zeros(ngnode,nper);

loadg_results=Array{Any}(nothing,nper+2,ngnode+1);
loadg_results[1:2,1]=["order" "gnode"];


for t in NPER
 
    loadg_results[t+2,1]=t;

    for n in NGN
        loadg_results[1:2,n+1]=[gnode[n, :order][1] gnode[n, :Node][1]];
     
        fl[n,t]=gc(n,t);

        loadg_results[t+2,n+1]=fl[n,t];              
    
     end
end




#Get r y x of branches
r(br)=branch[br,:r_pu];
x(br)=branch[br,:x_pu];
Capacity(br)=branch[br,:Pmax_MW];

#Get convex sub 0
b_min(pp)=0;
a_min(pp)=gnode[findall(gnode[:,:"Node"].==pipe[pp,:First])[1],:Pre_min]+gnode[findall(gnode[:,:"Node"].==pipe[pp,:Terminal])[1],:Pre_min]*pipe[pp,:sqrt_es];
b_max(pp)=gnode[findall(gnode[:,:"Node"].==pipe[pp,:First])[1],:Pre_max]-gnode[findall(gnode[:,:"Node"].==pipe[pp,:Terminal])[1],:Pre_min]*pipe[pp,:sqrt_es];
fg_max(pp)=sqrt((pipe[pp,:Cmn]/stageG)^2*(gnode[findall(gnode[:,:"Node"].==pipe[pp,:First])[1],:Pre_max]^2-(gnode[findall(gnode[:,:"Node"].==pipe[pp,:Terminal])[1],:Pre_min]*pipe[pp,:sqrt_es])^2));
fg_min(pp)=0;
a_max(pp)=gnode[findall(gnode[:,:"Node"].==pipe[pp,:First])[1],:Pre_max]+gnode[findall(gnode[:,:"Node"].==pipe[pp,:Terminal])[1],:Pre_max]*pipe[pp,:sqrt_es];

#Function that find the branch if (m,n) belongs to branches matrix
function map(m,n,branches)
    a=[m n]; 
    found=0;
    rowfound=0;
    Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
    for i in 1:size(Searchi,1)
        if Searchi[i,1] == 1
            found=1
            rowfound=i
        end
    end
    if found==0
        a=[n m];
        Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1) ]
        for i in 1:size(Searchi,1)
            if Searchi[i,1] == 1
                found=1
                rowfound=i
            end
        end
    end
    return rowfound
end

#Function that verify if the branch (m,n) exists
function exist(m,n,branches)
    a=[m n]; 
    found=0;
    rowfound=0;
    Searchi=Int[ a == [branches[i,1] branches[i,2]] for i=1:size(branches,1)]
    for i in 1:size(Searchi,1)
        if Searchi[i,1] == 1
            found=1
            rowfound=i
        end
    end
  
    return found
end

lmpp=zeros(nbus,nper);
#Power Flow Model

#PF=JuMP.Model(Juniper)
Fgmin=zeros(npipe,nper);
Fgmax=zeros(npipe,nper);
Amin=zeros(npipe,nper);
Amax=zeros(npipe,nper);
Bmin=zeros(npipe,nper);
Bmax=zeros(npipe,nper);

for pp in NPP, t in NPER
Amin[pp,t]=a_min(pp);
Fgmin[pp,t]=fg_min(pp);
Fgmax[pp,t]=fg_max(pp);
Amax[pp,t]=a_max(pp);
Bmax[pp,t]=b_max(pp);
Bmin[pp,t]=b_min(pp);
end

if con_ini=="ini"
    tolmip=50e-3
else
    tolmip=10e-3
end

function GP_CPk(Fgmin,Fgmax,Amin,Amax,Bmin,Bmax,k)
    println(string(""))
    println(string("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"))
    println(string("Iteration ",k," Model MISCOP"))

    PG = JuMP.Model(Gurobi.Optimizer);
    set_optimizer_attribute(PG, "MIPGap",tolmip)
    set_optimizer_attribute(PG, "Threads",thr)

    #Variables
    @variable(PG,Cost);
    #############Variables SP
    @variable(PG,0<=pt[NGEN,NPER]);# potencia termica
    @variable(PG,u[NGEN,NPER],binary=true);# binaria termica
    @variable(PG,0<=y[NGEN,NPER]<=1);# arranque termica
    @variable(PG,0<=w[NGEN,NPER]);# parada termica
    @variable(PG,0<=prer[NRER,NPER]);# potencia rer
    @variable(PG,0<=qthp[NGH,NPER]);# caudal turbinado hidro
    @variable(PG,uqthp[NGH,NPER],binary=true);# binaria hidro
    @variable(PG,0<=qvhp[NGH,NPER]);# caudal vertido hidro
    @variable(PG,0<=vol[NGHE,NPER]);# volumen hidro
    @variable(PG,f[NBR,NPER]);# flujo
    @variable(PG, -1.5708<=theta[NBUS, NPER]<=1.5708);# angulo electrico
    @variable(PG,0<=pshed[NBUS,NPER]);

    @variable(PG,0<=su[NWE,NPER]);
    #@variable(PG,cg[NGEN,NPER]);
    @variable(PG,0<=fg[NPP,NPER]);
    @variable(PG,fig[NPP,NPER]);
    @variable(PG,fog[NPP,NPER]);
    @variable(PG,0<=fco[NCO,NPER]);
    @variable(PG,0<=fre[NRE,NPER]);
    #@variable(PG,cco[NCO,NPER]);
    @variable(PG,lip[NPP,NPER]);
    @variable(PG,pr[NGN,NPER]);
    @variable(PG,0<=gshed[NGN,NPER]);
    #convex Variables
    #@variable(PG,a[NPP,NPER]);
    #@variable(PG,b[NPP,NPER]);
    @variable(PG,ka[NPP,NPER]);
    @variable(PG,la[NPP,NPER]);

    #CT
    @constraint(PG,p_max[g in NGEN, t in NPER],pt[g,t] <= plant[g,:Pmax]*u[g,t]*plant[g,:Disp]);#maximun power
    @constraint(PG,p_min[g in NGEN, t in NPER],pt[g,t] >= plant[g,:Pmin]*u[g,t]*plant[g,:Disp]);#minimun power
    @constraint(PG,ramp_g_max[g in NGEN, t in NPER],pt[g,t] - ( t==1 ? plant[g,:Pot_Ini] : pt[g,t-1])<=plant[g,:Ramp]);#Ramp max gen
    @constraint(PG,ramp_g_min[g in NGEN, t in NPER],pt[g,t] - ( t==1 ? plant[g,:Pot_Ini] : pt[g,t-1])>=-plant[g,:Ramp]);#Ramp min gen
    @constraint(PG, ecarr[g in NGEN, t in NPER], y[g,t] - w[g,t] == u[g,t] - ( t==1 ? plant[g,:Cond_I] : u[g,t-1]));# Arranque parada
    @constraint(PG, tminop[g in NGEN, tt in NPER], sum(y[g,t] for t in NPER if (t<=tt && (t>tt-plant[g,:TminO]))) <= u[g,tt] )# TMO
    @constraint(PG, tminfs[g in NGEN, tt in NPER], sum(w[g,t] for t in NPER if (t<=tt && (t>tt-plant[g,:TminF]))) <= 1-u[g,tt] )# TMFS
    @constraint(PG, ciclo_comb[cc in NCC, t in NPER], sum(ciclo[cc, :CCiclo]==plant[g, :CCiclo] ? u[g,t] : 0 for g in NGEN)<= 1 );# CC
    @constraint(PG, narr[g in NGEN, t in NPER], sum(y[g,t]) <= plant[g,:NA] )# Narranque
    @constraint(PG, nparr[g in NGEN, t in NPER], sum(w[g,t]) <= plant[g,:NP] )# Narranque
    @constraint(PG, forzado[g in NGEN, t in NPER], plant[g,:Forzado]<= u[g,t]);# Forzado
    @expression(PG,cg[g in NGEN,t in NPER],plant[g,:EffiD]*pt[g,t]+plant[g,:b]*u[g,t]);#Gas consum
    #@constraint(PG, tiniop[g in NGEN, t=1:(plant[g,:TminO]-plant[g,:T_ini_Op])], u[g,t]==1)
    #@constraint(PG, tinifs[g in NGEN, t=1:(plant[g,:TminF]-plant[g,:T_ini_Fs])], u[g,t]==0)
    #RER y RAC
    @constraint(PG,p_max_rer[gr in NRER, t in NPER],prer[gr,t] <= gerer[gr,t]);#maximun power rer
    @constraint(PG, psmax[n in NBUS, t in NPER],pshed[n,t] <= pcc[n,t]);#Rac power max
    #HIDRO
    @constraint(PG,php_max[gh in NGH, t in NPER],qthp[gh,t] <= planth[gh,:Pmax]*planth[gh,:Disp]/planth[gh,:Ren]*uqthp[gh,t]);#maximun q
    @constraint(PG,php_min[gh in NGH, t in NPER],qthp[gh,t] >= planth[gh,:Pmin]*planth[gh,:Disp]/planth[gh,:Ren]*uqthp[gh,t]);#minimun q
    @expression(PG,php[gh in NGH, t in NPER],planth[gh,:Ren]*qthp[gh,t]);#Redimiento
    @constraint(PG,toma_h[gh in NGHP, t in NPER],qthp[gh,t]+qvhp[gh,t] == planth[gh,:q]);#Bocatoma
    @constraint(PG,bal_vol[gh in NGHE, t in NPER],vol[gh,t]== ( t==1 ? planth[gh,:Vi] : vol[gh,t-1])
    +(planth[gh,:q]-qthp[gh,t]-qvhp[gh,t])*rel_h);#Balance Embalse
    @constraint(PG,min_vol[gh in NGHE, t in NPER],vol[gh,t]>= planth[gh,:Vmin])#Embalse min
    @constraint(PG,max_vol[gh in NGHE, t in NPER],vol[gh,t]<= planth[gh,:Vmax])#Embalse max
    @constraint(PG,goal_vol[gh in NGHE],vol[gh,nper]>= planth[gh,:Vf])#Embalse meta
    #TX
    @constraint(PG, reference[t in NPER], theta[1,t]==0);#Reference bus
    @constraint(PG, flowcap1[br in NBR, t in NPER],f[br,t]<= branch[br,:Pmax_MW]);#Capacity of the LT (in one direction)  
    @constraint(PG, flowcap2[br in NBR, t in NPER],f[br,t]>= -branch[br,:Pmax_MW]);#Capacity of the LT (in the opposite direction)
    @constraint(PG, pflow[br in NBR, t in NPER],
                        f[br,t]/100 ==  MVAbase*x(br)/(x(br)^2+r(br)^2)*(theta[mapbranch[br,1],t] - theta[mapbranch[br,2],t])/100 );#Flow   
    #Balance of power
    @constraint(PG, pbalance[n in NBUS, t in NPER], 
        deltaT*(sum( plant[g, :Node]==bus[n,:node] ? pt[g,t] : 0 for g in NGEN) 
        +sum( rer[gr, :node]==bus[n,:node] ? prer[gr,t] : 0 for gr in NRER )
        +sum( planth[gh, :node]==bus[n,:node] ? php[gh,t] : 0 for gh in NGH )
        + pshed[n,t] - pcc[n,t])
        == deltaT*sum( (exist(n,mapbranch[br,2],mapbranch[[br],:])==1) ? f[br,t] : 0 for br in NBR)  
           -deltaT*sum( (exist(mapbranch[br,1],n,mapbranch[[br],:])==1) ? f[br,t] : 0 for br in NBR));            

    ###############Constraint PG
    @constraint(PG,Cap_Well_min[w in NWE,t in NPER],su[w,t]>=well[w,:Sup_min_D]/stageG );#Capacity min well
    @constraint(PG,Cap_Well_max[w in NWE,t in NPER],su[w,t]<=well[w,:Sup_max_D]*well_cap[w,t]/stageG );#Capacity max well
    @constraint(PG, gsmax[gn in NGN, t in NPER],gshed[gn,t] <= fl[gn,t]);#Rac gas max
    @constraint(PG,pres_min[gn in NGN, t in NPER],pr[gn,t]>=gnode[gn,:Pre_min]);#Min Pression
    @constraint(PG,pres_max[gn in NGN, t in NPER],pr[gn,t]<=gnode[gn,:Pre_max]);#Max Pression
    ##COMPRESOR
    @constraint(PG,rel_min[com in NCO, t in NPER],pr[mapcom[com,1],t]*comp[com,:Rcomp_min]<=pr[mapcom[com,2],t]);#Rel min Pression
    @constraint(PG,rel_max[com in NCO, t in NPER],pr[mapcom[com,1],t]*comp[com,:Rcomp_max]>=pr[mapcom[com,2],t]);#Rel max Pression
    @expression(PG,cco[com in NCO, t in NPER],comp[com,:fcon]*(pr[mapcom[com,2],t]-pr[mapcom[com,1],t])/stageG);#Consumo de compresor
    ##REGULADOR DE PRESION
    @expression(PG,cre[re in NRE, t in NPER],fre[re,t]*reg[re,:fcon]);#Consumo de regulador de presion
    #@constraint(PG,out_reg_pres[re in NRE, t in NPER],pr[mapreg[re,2],t]==reg[re,:Pout]);#Presión de salida
    #@constraint(PG,int_reg_pmin[re in NRE, t in NPER],pr[mapreg[re,1],t]>=reg[re,:Pin_min]);#Presion mínima de entrada
    #@constraint(PG,int_reg_pmax[re in NRE, t in NPER],pr[mapreg[re,1],t]<=reg[re,:Pin_max]);#Presion mnima de entrada
    @constraint(PG,rel_minr[re in NRE, t in NPER],pr[mapreg[re,1],t]*reg[re,:Rreg_min]<=pr[mapreg[re,2],t]);#Rel min Pression
    @constraint(PG,rel_maxr[re in NRE, t in NPER],pr[mapreg[re,1],t]*reg[re,:Rreg_max]>=pr[mapreg[re,2],t]);#Rel max Pression
    #Balance of gas
    @constraint(PG,gbalance[gn in NGN,t in NPER],
                sum(well[w,:Gnode]==gnode[gn,:Node] ? su[w,t] : 0 for w in NWE)
                + gshed[gn,t]
                -sum(plant[g,:Gnode]==gnode[gn,:Node] ? cg[g,t] : 0 for g in NGEN)
                -fl[gn,t]
                -sum(comp[com,:First]==gnode[gn,:Node] ? cco[com,t] : 0 for com in NCO)
                -sum(reg[re,:Terminal]==gnode[gn,:Node] ? cre[re,t] : 0 for re in NRE)
                == sum((exist(gn,mappipe[pp,2],mappipe[[pp],:])==1) ? fig[pp,t] : 0 for pp in NPP)
                -sum((exist(mappipe[pp,1],gn,mappipe[[pp],:])==1) ? fog[pp,t] : 0 for pp in NPP)
                +sum((exist(gn,mapcom[com,2],mapcom[[com],:])==1) ? fco[com,t] : 0 for com in NCO)
                -sum((exist(mapcom[com,1],gn,mapcom[[com],:])==1) ? fco[com,t] : 0 for com in NCO)
                +sum((exist(gn,mapreg[re,2],mapreg[[re],:])==1) ? fre[re,t] : 0 for re in NRE)
                -sum((exist(mapreg[re,1],gn,mapreg[[re],:])==1) ? fre[re,t] : 0 for re in NRE));
    #Gas Flow Convexification|
    @constraint(PG, gflow[pp in NPP, t in NPER],
                        [ (pr[mappipe[pp,1],t]) ; [ (fg[pp,t]/(pipe[pp,:Cmn]/stageG)) , (pipe[pp,:sqrt_es]*pr[mappipe[pp,2],t]) , 0  ] ] in SecondOrderCone() );
    #=     @constraint(PG, gflow[pp in NPP, t in NPER],
                    (fg[pp,t]/(pipe[pp,:Cmn]/stageG))^2 <=  (pr[mappipe[pp,1],t])^2 - ((pipe[pp,:sqrt_es]*pr[mappipe[pp,2],t])^2 ); =#
    #= @constraint(PG, eqa[pp in NPP, t in NPER], a[pp,t] ==  pr[mappipe[pp,1],t] + pr[mappipe[pp,2],t]); =#
    @expression(PG, a[pp in NPP, t in NPER], pr[mappipe[pp,1],t] + pipe[pp,:sqrt_es]*pr[mappipe[pp,2],t]);
    #= @constraint(PG, eqb[pp in NPP, t in NPER], b[pp,t] ==  pr[mappipe[pp,1],t] - pr[mappipe[pp,2],t]); =#
    @expression(PG, b[pp in NPP, t in NPER], pr[mappipe[pp,1],t] - pipe[pp,:sqrt_es]*pr[mappipe[pp,2],t]);
    @constraint(PG, eqI[pp in NPP, t in NPER], 100*(ka[pp,t]*(1/(pipe[pp,:Cmn]/stageG))^2) >= 100*(la[pp,t])); 
    #@constraint(PG, eqII[pp in NPP, t in NPER],  ka[pp,t] >= fg[pp,t]^2 );
    @constraint(PG, eqII[pp in NPP, t in NPER], [ ka[pp,t]+1 ; [ 2*fg[pp,t] , ka[pp,t]-1 ] ] in SecondOrderCone());
    @constraint(PG, eq1[pp in NPP, t in NPER],(ka[pp,t]) <= ((Fgmax[pp,t]+Fgmin[pp,t])*fg[pp,t]-Fgmax[pp,t]*Fgmin[pp,t]));
    @constraint(PG, eq2[pp in NPP, t in NPER],la[pp,t] >= Amin[pp,t]*b[pp,t]+Bmin[pp,t]*a[pp,t]-Amin[pp,t]*Bmin[pp,t]);
    @constraint(PG, eq3[pp in NPP, t in NPER],(la[pp,t])/100>= (Amax[pp,t]*b[pp,t]+Bmax[pp,t]*a[pp,t]-Amax[pp,t]*Bmax[pp,t])/100);
    @constraint(PG, eq4[pp in NPP, t in NPER],(la[pp,t])/100 <= (Amin[pp,t]*b[pp,t]+Bmax[pp,t]*a[pp,t]-Amin[pp,t]*Bmax[pp,t])/100);
    @constraint(PG, eq5[pp in NPP, t in NPER],la[pp,t] <= (Amax[pp,t]*b[pp,t]+Bmin[pp,t]*a[pp,t]-Amax[pp,t]*Bmin[pp,t]));
    #Gas in out Flow
    @constraint(PG, g_io_flow[pp in NPP, t in NPER],
                        fg[pp,t]== (fig[pp,t] + fog[pp,t])*0.5 );
    #Gas line pack and flow
    if con_ini=="ini"
        @constraint(PG, lp_flow[pp in NPP, t in NPER],0 == fig[pp,t] - fog[pp,t]);  
        #@constraint(PG, lp_ini, 28.3168*lpini == sum(pipe[pp,:Linepack]==1 ? lip[pp,1] : 0 for pp in NPP)); #nesw for CI
    else
        #@constraint(PG, nlp_flow[pp in NPPNL, t in NPER],0 == fig[pp,t] - fog[pp,t]);  
        @constraint(PG, lp_flow[pp in NPP, t in NPER],lip[pp,t]-( t==1 ? pipe[pp,:L_0] : lip[pp,t-1] ) == fig[pp,t] - fog[pp,t]);    
        @constraint(PG, lp_goal[pp in NPP],lip[pp,nper]>= pipe[pp,:L_f] );#Gas line pack meta
        #@constraint(PG, lp_goal2[pp in NPP],lip[pp,nper]<= 1.02*pipe[pp,:L_f] );#Gas line pack meta
    end               
    @constraint(PG, linepack[pp in NPP, t in NPER],lip[pp,t]==(pipe[pp,:Kmn])/2*(pr[mappipe[pp,1],t]+pr[mappipe[pp,2],t]));#Line Pack
    #Gas in out Flow
    @constraint(PG, limite_flow[pp in NPP, t in NPER],fg[pp,t]<= pipe[pp,:Cap]/stageG);
    # Operación Planta_no varia
    #@constraint(PG, op_well_1[t in NPER],su[1,t]/10000<= (9.9249*pr[1,t]-77216)/stageG/10000);
    #@constraint(PG, op_well_2[t in NPER],su[1,t]/10000<= (4.0072*pr[1,t]-12535)/stageG/10000);
    #@constraint(PG, op_well_3[t in NPER],su[1,t]/10000>= (29.25*pr[1,t]-391175)/stageG/10000);
    #@constraint(PG, op_well_4[t in NPER],su[1,t]/10000>= (9.4493*pr[1,t]-104068)/stageG/10000);
    #@constraint(PG, op_well_5[t in NPER],su[1,t]/10000>= (2.4035*pr[1,t]-11027)/stageG/10000);

    #Objective
    @objective(PG, Min, Cost);
    #Cost           
    if con_ini=="ini"
#=         @constraint(PG,Cost==(deltaT*sum(plant[g,:cvnc]*pt[g,t] for g in NGEN,t in NPER)+sum(plant[g,:Co]*cg[g,t] for g in NGEN,t in NPER)
        +sum(plant[g,:C_start]*y[g,t] for g in NGEN,t in NPER)+sum(plant[g,:C_down]*w[g,t] for g in NGEN,t in NPER)
        +deltaT*sum(pshed[n,t]*price_shed for n in NBUS, t in NPER)
        +deltaT*sum(planth[gh,:cost]*php[gh,t] for gh in NGH,t in NPER)
        +sum(well[w,:C_sup_D]*su[w,t] for w in NWE,t in NPER)
        +sum(gshed[gn,t]*g_price_shed for gn in NGN, t in NPER))); =#
        @constraint(PG,Cost==(deltaT*sum(plant[g,:cvnc]*pt[g,t] for g in NGEN,t in NPER)+sum(plant[g,:Co]*cg[g,t] for g in NGEN,t in NPER)
        +sum(plant[g,:C_start]*y[g,t] for g in NGEN,t in NPER)+sum(plant[g,:C_down]*w[g,t] for g in NGEN,t in NPER)
        +deltaT*sum(pshed[n,t]*price_shed for n in NBUS, t in NPER)
        +deltaT*sum(planth[gh,:cost]*php[gh,t] for gh in NGH,t in NPER)
        +sum(well[w,:C_sup_D]*su[w,t] for w in NWE,t in NPER)
        +sum(gshed[gn,t]*g_price_shed for gn in NGN, t in NPER)-sum(lip[pp,1]*11 for pp in NPP)));
    else
        @constraint(PG,Cost==(deltaT*sum(plant[g,:cvnc]*pt[g,t] for g in NGEN,t in NPER)+sum(plant[g,:Co]*cg[g,t] for g in NGEN,t in NPER)
        +sum(plant[g,:C_start]*y[g,t] for g in NGEN,t in NPER)+sum(plant[g,:C_down]*w[g,t] for g in NGEN,t in NPER)
        +deltaT*sum(pshed[n,t]*price_shed for n in NBUS, t in NPER)
        +deltaT*sum(planth[gh,:cost]*php[gh,t] for gh in NGH,t in NPER)
        +sum(well[w,:C_sup_D]*su[w,t] for w in NWE,t in NPER)
        +sum(gshed[gn,t]*g_price_shed for gn in NGN, t in NPER)));
    end   


    optimize!(PG);

    return  JuMP.value.(fg),        #1
            JuMP.value.(a),         #2
            JuMP.value.(b),         #3
            JuMP.value.(pr),        #4
            JuMP.value.(Cost),      #5
            JuMP.value.(pt),        #6
            JuMP.value.(su),        #7
            JuMP.value.(cg),        #8
            JuMP.value.(y),         #9
            JuMP.value.(w),         #10
            JuMP.value.(gshed),     #11
            JuMP.value.(pshed),     #12
            JuMP.value.(f),         #13
            JuMP.value.(lip),       #14
            JuMP.value.(fig),       #15
            JuMP.value.(fog),       #16 
            JuMP.value.(cco),       #17
            JuMP.value.(fco),       #18
            JuMP.value.(prer),      #19
            JuMP.value.(php),       #20
            JuMP.value.(vol),       #21
            JuMP.value.(fre),       #22
            JuMP.value.(cre),       #23
			JuMP.value.(u),       	#24
            JuMP.value.(uqthp)      #25
			
    end

#ee=[0.5 0.35 0.4 0.25 0.125 0.075];
ee=[0.5]
iter=size(ee,2);
global k
k=1;

global error_mean=ones(iter)
global costs=Array{Any}(nothing,iter);
global error_max=Array{Any}(nothing,iter);
start_time=Base.time_ns();
global solution=GP_CPk(Fgmin,Fgmax,Amin,Amax,Bmin,Bmax,k); 

    # Array of violations
    global viol=Array{Any}(nothing,npipe,nper);
    global violabs=Array{Any}(nothing,npipe,nper);
    for pp in NPP, t in NPER
        viol[pp,t]=[((solution[4][mappipe[pp,1],t])^2 - (pipe[pp,:sqrt_es]*solution[4][mappipe[pp,2],t])^2-(solution[1][pp,t])^2/(pipe[pp,:Cmn]/stageG)^2)/((solution[4][mappipe[pp,1],t])^2)];
        violabs[pp,t]=[abs((solution[4][mappipe[pp,1],t])^2 - (pipe[pp,:sqrt_es]*solution[4][mappipe[pp,2],t])^2-(solution[1][pp,t])^2/(pipe[pp,:Cmn]/stageG)^2)/((solution[4][mappipe[pp,1],t])^2)];
    end

    global ViolmaxT=maximum(maximum(violabs,dims=1));  
    global Violmean1=sum(viol)/nper/npipe;
    error_mean[k]=Violmean1[1];
    error_max[k]=ViolmaxT[1,1];
    costs[k]=solution[5];
    global fmink=(1-ee[k])*(solution[1]);
    global fmaxk=(1+ee[k])*(solution[1]);
    global amink=(1-ee[k])*(solution[2]);
    global amaxk=(1+ee[k])*(solution[2]);
    global bmink=(1-ee[k])*(solution[3]);
    global bmaxk=(1+ee[k])*(solution[3]);
global k
k=2;

while k<=iter && abs(error_mean[k-1])>=0.02

        #while abs(error_mean[k-1])>=0.05
        #solution=GP_CPk(fmink,fmaxk,amink,amaxk,bmink,bmaxk);
    global    solution=GP_CPk(Fgmin,fmaxk,amink,Amax,Bmin,bmaxk,k); 
        #solution=GP_CPk(FF_min,FF_max,AA_min,amaxk,BB_min,bmaxk);

    global fmink=(1-ee[k])*(solution[1]);
    global fmaxk=(1+ee[k])*(solution[1]);
    global amink=(1-ee[k])*(solution[2]);
    global amaxk=(1+ee[k])*(solution[2]);
    global bmink=(1-ee[k])*(solution[3]);
    global bmaxk=(1+ee[k])*(solution[3]);
    
    # Array of violations
    global viol=Array{Any}(nothing,npipe,nper);
    global violabs=Array{Any}(nothing,npipe,nper);
    for pp in NPP, t in NPER
        viol[pp,t]=[((solution[4][mappipe[pp,1],t])^2 - (pipe[pp,:sqrt_es]*solution[4][mappipe[pp,2],t])^2-(solution[1][pp,t])^2/(pipe[pp,:Cmn]/stageG)^2)/((solution[4][mappipe[pp,1],t])^2)];
        violabs[pp,t]=[abs((solution[4][mappipe[pp,1],t])^2 - (pipe[pp,:sqrt_es]*solution[4][mappipe[pp,2],t])^2-(solution[1][pp,t])^2/(pipe[pp,:Cmn]/stageG)^2)/((solution[4][mappipe[pp,1],t])^2)];
    end

    global ViolmaxT=maximum(maximum(violabs,dims=1));  
    global Violmean1=sum(viol)/nper/npipe;
    error_mean[k]=Violmean1[1];
    error_max[k]=ViolmaxT[1,1];
    costs[k]=solution[5];

    global k+=1
end

    #solutions
    global fg_cx=solution[1]
    global a_cx=solution[2]
    global b_cx=solution[3]
    global pr_cx=solution[4]
    global Cost_cx=solution[5]
    global pt_cx=solution[6]
    global su_cx=solution[7]
    global cg_cx=solution[8]
    global y_cx=solution[9]
    global w_cx=solution[10]
    global gshed_cx=solution[11]
    global pshed_cx=solution[12]
    global f_cx=solution[13]
    global lip_cx=solution[14]
    global fig_cx=solution[15]
    global fog_cx=solution[16]
    global cco_cx=solution[17]
    global fco_cx=solution[18]
    global prer_cx=solution[19]
    global php_cx=solution[20]
    global vol_cx=solution[21]
    global fre_cx=solution[22]
    global cre_cx=solution[23]
    global u_cx=solution[24]
    global uqthp_cx=solution[25]

iter=k-1

println(string(""))
println(string("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"))
println(string("Iteration ModelCMG  MISCOP"))

CPG = JuMP.Model(Gurobi.Optimizer);
set_optimizer_attribute(CPG, "Threads",thr)
set_optimizer_attribute(CPG, "MIPGap",6.0e-3)

#Variables
@variable(CPG,Cost);
#############Variables SP
@variable(CPG,0<=pt[NGEN,NPER]);# potencia termica
@variable(CPG,0<=prer[NRER,NPER]);# potencia rer
@variable(CPG,0<=qthp[NGH,NPER]);# caudal turbinado hidro
@variable(CPG,0<=qvhp[NGH,NPER]);# caudal vertido hidro
@variable(CPG,0<=vol[NGHE,NPER]);# volumen hidro
@variable(CPG,f[NBR,NPER]);# flujo
@variable(CPG, -1.5708<=theta[NBUS, NPER]<=1.5708);# angulo electrico

@variable(CPG,0<=su[NWE,NPER]);
#@variable(CPG,cg[NGEN,NPER]);
@variable(CPG,0<=fg[NPP,NPER]);
@variable(CPG,fig[NPP,NPER]);
@variable(CPG,fog[NPP,NPER]);
@variable(CPG,0<=fco[NCO,NPER]);
@variable(CPG,0<=fre[NRE,NPER]);
#@variable(CPG,cco[NCO,NPER]);
@variable(CPG,lip[NPP,NPER]);
@variable(CPG,pr[NGN,NPER]);
@variable(CPG,0<=pshed[NBUS,NPER]);
@variable(CPG,0<=gshed[NGN,NPER]);
#convex Variables
@variable(CPG,a[NPP,NPER]);
@variable(CPG,b[NPP,NPER]);
@variable(CPG,ka[NPP,NPER]);
@variable(CPG,la[NPP,NPER]);

#CT
@expression(CPG, u[g in NGEN, t in NPER], u_cx[g,t]);
@expression(CPG, w[g in NGEN, t in NPER], w_cx[g,t]);
@expression(CPG, y[g in NGEN, t in NPER], y_cx[g,t]);
@expression(CPG,uqthp[gh in NGH, t in NPER], uqthp_cx[gh,t]);
@constraint(CPG,p_max[g in NGEN, t in NPER],pt[g,t] <= plant[g,:Pmax]*u[g,t]*plant[g,:Disp]);#maximun power
@constraint(CPG,p_min[g in NGEN, t in NPER],pt[g,t] >= plant[g,:Pmin]*u[g,t]*plant[g,:Disp]);#minimun power
@constraint(CPG,ramp_g_max[g in NGEN, t in NPER],pt[g,t] - ( t==1 ? plant[g,:Pot_Ini] : pt[g,t-1])<=plant[g,:Ramp]);#Ramp max gen
@constraint(CPG,ramp_g_min[g in NGEN, t in NPER],pt[g,t] - ( t==1 ? plant[g,:Pot_Ini] : pt[g,t-1])>=-plant[g,:Ramp]);#Ramp min gen
@constraint(CPG, ciclo_comb[cc in NCC, t in NPER], sum(ciclo[cc, :CCiclo]==plant[g, :CCiclo] ? u[g,t] : 0 for g in NGEN)<= 1 );# CC
@expression(CPG,cg[g in NGEN,t in NPER],plant[g,:EffiD]*pt[g,t]+plant[g,:b]*u[g,t]);#Gas consum
#@constraint(CPF, tiniop[g in NGEN, t=1:(plant[g,:TminO]-plant[g,:T_ini_Op])], u[g,t]==1)
#@constraint(CPF, tinifs[g in NGEN, t=1:(plant[g,:TminF]-plant[g,:T_ini_Fs])], u[g,t]==0)
#RER y RAC
@constraint(CPG,p_max_rer[gr in NRER, t in NPER],prer[gr,t] <= gerer[gr,t]);#maximun power rer
@constraint(CPG, psmax[n in NBUS, t in NPER],pshed[n,t] <= pcc[n,t]);#Rac power max
#HIDRO
@constraint(CPG,php_max[gh in NGH, t in NPER],qthp[gh,t] <= planth[gh,:Pmax]*planth[gh,:Disp]/planth[gh,:Ren]*uqthp[gh,t]);#maximun q
@constraint(CPG,php_min[gh in NGH, t in NPER],qthp[gh,t] >= planth[gh,:Pmin]*planth[gh,:Disp]/planth[gh,:Ren]*uqthp[gh,t]);#minimun q
@expression(CPG,php[gh in NGH, t in NPER],planth[gh,:Ren]*qthp[gh,t]);#Redimiento
@constraint(CPG,toma_h[gh in NGHP, t in NPER],qthp[gh,t]+qvhp[gh,t] == planth[gh,:q]);#Bocatoma
@constraint(CPG,bal_vol[gh in NGHE, t in NPER],vol[gh,t]== ( t==1 ? planth[gh,:Vi] : vol[gh,t-1])
    +(planth[gh,:q]-qthp[gh,t]-qvhp[gh,t])*rel_h);#Balance Embalse
@constraint(CPG,min_vol[gh in NGHE, t in NPER],vol[gh,t]>= planth[gh,:Vmin])#Embalse min
@constraint(CPG,max_vol[gh in NGHE, t in NPER],vol[gh,t]<= planth[gh,:Vmax])#Embalse max
@constraint(CPG,goal_vol[gh in NGHE],vol[gh,nper]>= planth[gh,:Vf])#Embalse meta
#TX
@constraint(CPG, reference[t in NPER], theta[1,t]==0);#Reference bus
@constraint(CPG, flowcap1[br in NBR, t in NPER],f[br,t]<= branch[br,:Pmax_MW]);#Capacity of the LT (in one direction)  
@constraint(CPG, flowcap2[br in NBR, t in NPER],f[br,t]>= -branch[br,:Pmax_MW]);#Capacity of the LT (in the opposite direction)
@constraint(CPG, CPGlow[br in NBR, t in NPER],
                        f[br,t]/100 ==  MVAbase*x(br)/(x(br)^2+r(br)^2)*(theta[mapbranch[br,1],t] - theta[mapbranch[br,2],t])/100 );#Flow   
#Balance of power
@constraint(CPG, pbalance[n in NBUS, t in NPER], 
        deltaT*(sum( plant[g, :Node]==bus[n,:node] ? pt[g,t] : 0 for g in NGEN) 
        +sum( rer[gr, :node]==bus[n,:node] ? prer[gr,t] : 0 for gr in NRER )
        +sum( planth[gh, :node]==bus[n,:node] ? php[gh,t] : 0 for gh in NGH )
        + pshed[n,t] - pcc[n,t])
        == deltaT*sum( (exist(n,mapbranch[br,2],mapbranch[[br],:])==1) ? f[br,t] : 0 for br in NBR)  
           -deltaT*sum( (exist(mapbranch[br,1],n,mapbranch[[br],:])==1) ? f[br,t] : 0 for br in NBR));            

###############Constraint CPG
@constraint(CPG,Cap_Well_min[w in NWE,t in NPER],su[w,t]>=well[w,:Sup_min_D]/stageG );#Capacity min well
@constraint(CPG,Cap_Well_max[w in NWE,t in NPER],su[w,t]<=well[w,:Sup_max_D]*well_cap[w,t]/stageG );#Capacity max well
@constraint(CPG, gsmax[gn in NGN, t in NPER],gshed[gn,t] <= fl[gn,t]);#Rac gas max
@constraint(CPG,pres_min[gn in NGN, t in NPER],pr[gn,t]>=gnode[gn,:Pre_min]);#Min Pression
@constraint(CPG,pres_max[gn in NGN, t in NPER],pr[gn,t]<=gnode[gn,:Pre_max]);#Max Pression
##COMPRESOR
@constraint(CPG,rel_min[com in NCO, t in NPER],pr[mapcom[com,1],t]*comp[com,:Rcomp_min]<=pr[mapcom[com,2],t]);#Rel min Pression
@constraint(CPG,rel_max[com in NCO, t in NPER],pr[mapcom[com,1],t]*comp[com,:Rcomp_max]>=pr[mapcom[com,2],t]);#Rel max Pression
@expression(CPG,cco[com in NCO, t in NPER],comp[com,:fcon]*(pr[mapcom[com,2],t]-pr[mapcom[com,1],t])/stageG);#Consumo de compresor
##REGULADOR DE PRESION
@expression(CPG,cre[re in NRE, t in NPER],fre[re,t]*reg[re,:fcon]);#Consumo de regulador de presion
#@constraint(CPG,out_reg_pres[re in NRE, t in NPER],pr[mapreg[re,2],t]==reg[re,:Pout]);#Presión de salida
#@constraint(CPG,int_reg_pmin[re in NRE, t in NPER],pr[mapreg[re,1],t]>=reg[re,:Pin_min]);#Presion mínima de entrada
#@constraint(CPG,int_reg_pmax[re in NRE, t in NPER],pr[mapreg[re,1],t]<=reg[re,:Pin_max]);#Presion mnima de entrada
@constraint(CPG,rel_minr[re in NRE, t in NPER],pr[mapreg[re,1],t]*reg[re,:Rreg_min]<=pr[mapreg[re,2],t]);#Rel min Pression
@constraint(CPG,rel_maxr[re in NRE, t in NPER],pr[mapreg[re,1],t]*reg[re,:Rreg_max]>=pr[mapreg[re,2],t]);#Rel max Pression
#Balance of gas
@constraint(CPG,gbalance[gn in NGN,t in NPER],
            (sum(well[w,:Gnode]==gnode[gn,:Node] ? su[w,t] : 0 for w in NWE)
            + gshed[gn,t]
            -sum(plant[g,:Gnode]==gnode[gn,:Node] ? cg[g,t] : 0 for g in NGEN)
            -fl[gn,t]
            -sum(comp[com,:First]==gnode[gn,:Node] ? cco[com,t] : 0 for com in NCO)
            -sum(reg[re,:Terminal]==gnode[gn,:Node] ? cre[re,t] : 0 for re in NRE))*1000
            == (sum((exist(gn,mappipe[pp,2],mappipe[[pp],:])==1) ? fig[pp,t] : 0 for pp in NPP)
                -sum((exist(mappipe[pp,1],gn,mappipe[[pp],:])==1) ? fog[pp,t] : 0 for pp in NPP)
                +sum((exist(gn,mapcom[com,2],mapcom[[com],:])==1) ? fco[com,t] : 0 for com in NCO)
                -sum((exist(mapcom[com,1],gn,mapcom[[com],:])==1) ? fco[com,t] : 0 for com in NCO)
                +sum((exist(gn,mapreg[re,2],mapreg[[re],:])==1) ? fre[re,t] : 0 for re in NRE)
                -sum((exist(mapreg[re,1],gn,mapreg[[re],:])==1) ? fre[re,t] : 0 for re in NRE))*1000);

#Gas Flow Linealization
@constraint(CPG, gflow[pp in NPP, t in NPER],
                    (fg_cx[pp,t]^2+2*fg_cx[pp,t]*(fg[pp,t]-fg_cx[pp,t]))/(pipe[pp,:Cmn]/stageG)^2 ==  
					(pr_cx[mappipe[pp,1],t])^2 +2*pr_cx[mappipe[pp,1],t]*(pr[mappipe[pp,1],t]-pr_cx[mappipe[pp,1],t])
					-((pipe[pp,:sqrt_es]*pr_cx[mappipe[pp,2],t])^2+ 2*pipe[pp,:sqrt_es]*pr_cx[mappipe[pp,2],t]*(pipe[pp,:sqrt_es]*pr[mappipe[pp,2],t]-pipe[pp,:sqrt_es]*pr_cx[mappipe[pp,2],t])));
#Gas in out Flow

@constraint(CPG, g_io_flow[pp in NPP, t in NPER],
                        fg[pp,t]== (fig[pp,t] + fog[pp,t])*0.5 );
    if con_ini=="ini"
        @constraint(CPG, lp_flow[pp in NPP, t in NPER],0 == fig[pp,t] - fog[pp,t]);  
        #@constraint(CPG, lp_ini, 28.3168*lpini == sum(pipe[pp,:Linepack]==1 ? lip[pp,1] : 0 for pp in NPP)); #nesw for CI
    else
        #@constraint(CPG, nlp_flow[pp in NPPNL, t in NPER],0 == fig[pp,t] - fog[pp,t]);  
        @constraint(CPG, lp_flow[pp in NPP, t in NPER],lip[pp,t]-( t==1 ? pipe[pp,:L_0] : lip[pp,t-1] ) == fig[pp,t] - fog[pp,t]);    
        @constraint(CPG, lp_goal[pp in NPP],lip[pp,nper]>= pipe[pp,:L_f] );#Gas line pack meta
        #@constraint(CPG, lp_goal2[pp in NPP],lip[pp,nper]<= 1.02*pipe[pp,:L_f] );#Gas line pack meta
    end                  
@constraint(CPG, linepack[pp in NPP, t in NPER],lip[pp,t]==(pipe[pp,:Kmn])/2*(pr[mappipe[pp,1],t]+pr[mappipe[pp,2],t]));#Line Pack
#Gas in out Flow
#Gas in out Flow
@constraint(CPG, limite_flow[pp in NPP, t in NPER],fg[pp,t]<= pipe[pp,:Cap]/stageG);
# Operación Planta_no varia
#@constraint(CPG, op_well_1[t in NPER],su[1,t]/10000<= (9.9249*pr[1,t]-77216)/stageG/10000);
#@constraint(CPG, op_well_2[t in NPER],su[1,t]/10000<= (4.0072*pr[1,t]-12535)/stageG/10000);
#@constraint(CPG, op_well_3[t in NPER],su[1,t]/10000>= (29.25*pr[1,t]-391175)/stageG/10000);
#@constraint(CPG, op_well_4[t in NPER],su[1,t]/10000>= (9.4493*pr[1,t]-104068)/stageG/10000);
#@constraint(CPG, op_well_5[t in NPER],su[1,t]/10000>= (2.4035*pr[1,t]-11027)/stageG/10000);

#Objective
@objective(CPG, Min, Cost);
#Cost           
@constraint(CPG,Cost==(deltaT*sum(plant[g,:cvnc]*pt[g,t] for g in NGEN,t in NPER)+sum(plant[g,:Co]*cg[g,t] for g in NGEN,t in NPER)
                    +sum(plant[g,:C_start]*y[g,t] for g in NGEN,t in NPER)+sum(plant[g,:C_down]*w[g,t] for g in NGEN,t in NPER)
                    +deltaT*sum(pshed[n,t]*price_shed for n in NBUS, t in NPER)
                    +deltaT*sum(planth[gh,:cost]*php[gh,t] for gh in NGH,t in NPER)
                    +sum(well[w,:C_sup_D]*su[w,t] for w in NWE,t in NPER)
                    +sum(gshed[gn,t]*g_price_shed for gn in NGN, t in NPER)));


optimize!(CPG);

global lmpp_a=JuMP.dual.(pbalance);
global lfl1_a=JuMP.dual.(flowcap1);
global lfl2_a=JuMP.dual.(flowcap2);
global lmpp_gn=JuMP.dual.(gbalance);
global lp_min_gn=JuMP.dual.(pres_min);
global lp_max_gn=JuMP.dual.(pres_max); 
global lp_flow_gn=JuMP.dual.(limite_flow); 

#Iteración para Cmg solo SP
bolsagn=sum(sum(plant[g,:Bolsa]=="B_Camisea" ? cg_cx[g,t] : 0 for g in NGEN) for t in NPER)
matriz = reshape(cg_cx, (60, nper));
sum_col=sum(matriz, dims=2);
global cg_bolsa_d=vec(sum_col);
ND=1:ndays;
println(string(""))
println(string("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"))
println(string("Iteration Model Cmg SP"))
CPF = JuMP.Model(Gurobi.Optimizer);
set_optimizer_attribute(CPF, "MIPGap",5.00e-3)
set_optimizer_attribute(CPF, "Threads",thr)
#############Variables SP
@variable(CPF,0<=pt[NGEN,NPER]);# potencia CT
@variable(CPF,CostE);
#RER y RAC
@variable(CPF,0<=prer[NRER,NPER]);# potencia RER
@variable(CPF,0<=pshed[NBUS,NPER]);# racionamiento SP
#Tx
@variable(CPF, -1.5708<=theta[NBUS, NPER]<=1.5708);# angulo electrico
@variable(CPF,f[NBR,NPER]);# flujo de potencia
#CH
@variable(CPF,0<=qthp[NGH,NPER]);
@variable(CPF,0<=qvhp[NGH,NPER]);
@variable(CPF,0<=vol[NGHE,NPER]);
#############Restricciones SP
#CT
@expression(CPF, u[g in NGEN, t in NPER], u_cx[g,t]);
@expression(CPF, w[g in NGEN, t in NPER], w_cx[g,t]);
@expression(CPF, y[g in NGEN, t in NPER], y_cx[g,t]);
@expression(CPF,uqthp[gh in NGH, t in NPER], uqthp_cx[gh,t]);
@constraint(CPF,p_max[g in NGEN, t in NPER],pt[g,t] <= plant[g,:Pmax]*u[g,t]*plant[g,:Disp]);#maximun power
@constraint(CPF,p_min[g in NGEN, t in NPER],pt[g,t] >= plant[g,:Pmin]*u[g,t]*plant[g,:Disp]);#minimun power
@constraint(CPF,ramp_g_max[g in NGEN, t in NPER],pt[g,t] - ( t==1 ? plant[g,:Pot_Ini] : pt[g,t-1])<=plant[g,:Ramp]);#Ramp max gen
@constraint(CPF,ramp_g_min[g in NGEN, t in NPER],pt[g,t] - ( t==1 ? plant[g,:Pot_Ini] : pt[g,t-1])>=-plant[g,:Ramp]);#Ramp min gen
@constraint(CPF, ciclo_comb[cc in NCC, t in NPER], sum(ciclo[cc, :CCiclo]==plant[g, :CCiclo] ? u[g,t] : 0 for g in NGEN)<= 1 );# CC
@expression(CPF,cg[g in NGEN,t in NPER],plant[g,:EffiD]*pt[g,t]+plant[g,:b]*u[g,t]);#Gas consum
#@constraint(CPF, tiniop[g in NGEN, t=1:(plant[g,:TminO]-plant[g,:T_ini_Op])], u[g,t]==1)
#@constraint(CPF, tinifs[g in NGEN, t=1:(plant[g,:TminF]-plant[g,:T_ini_Fs])], u[g,t]==0)
#RER y RAC
@constraint(CPF,p_max_rer[gr in NRER, t in NPER],prer[gr,t] <= gerer[gr,t]);#maximun power rer
@constraint(CPF, psmax[n in NBUS, t in NPER],pshed[n,t] <= pcc[n,t]);#Rac power max
#HIDRO
@constraint(CPF,php_max[gh in NGH, t in NPER],qthp[gh,t] <= planth[gh,:Pmax]*planth[gh,:Disp]/planth[gh,:Ren]*uqthp[gh,t]);#maximun q
@constraint(CPF,php_min[gh in NGH, t in NPER],qthp[gh,t] >= planth[gh,:Pmin]*planth[gh,:Disp]/planth[gh,:Ren]*uqthp[gh,t]);#minimun q
@expression(CPF,php[gh in NGH, t in NPER],planth[gh,:Ren]*qthp[gh,t]);#Redimiento
@constraint(CPF,toma_h[gh in NGHP, t in NPER],qthp[gh,t]+qvhp[gh,t] == planth[gh,:q]);#Bocatoma
@constraint(CPF,bal_vol[gh in NGHE, t in NPER],vol[gh,t]== ( t==1 ? planth[gh,:Vi] : vol[gh,t-1])
    +(planth[gh,:q]-qthp[gh,t]-qvhp[gh,t])*rel_h);#Balance Embalse
@constraint(CPF,min_vol[gh in NGHE, t in NPER],vol[gh,t]>= planth[gh,:Vmin])#Embalse min
@constraint(CPF,max_vol[gh in NGHE, t in NPER],vol[gh,t]<= planth[gh,:Vmax])#Embalse max
@constraint(CPF,goal_vol[gh in NGHE],vol[gh,nper]>= planth[gh,:Vf])#Embalse meta
#TX
@constraint(CPF, reference[t in NPER], theta[1,t]==0);#Reference bus
@constraint(CPF, flowcap1[br in NBR, t in NPER],f[br,t]<= branch[br,:Pmax_MW]);#Capacity of the LT (in one direction)  
@constraint(CPF, flowcap2[br in NBR, t in NPER],f[br,t]>= -branch[br,:Pmax_MW]);#Capacity of the LT (in the opposite direction)
@constraint(CPF, CPFlow[br in NBR, t in NPER],
                        f[br,t]/100 ==  MVAbase*x(br)/(x(br)^2+r(br)^2)*(theta[mapbranch[br,1],t] - theta[mapbranch[br,2],t])/100 );#Flow   
#Balance of power
@constraint(CPF, pbalance[n in NBUS, t in NPER], 
        deltaT*(sum( plant[g, :Node]==bus[n,:node] ? pt[g,t] : 0 for g in NGEN) 
        +sum( rer[gr, :node]==bus[n,:node] ? prer[gr,t] : 0 for gr in NRER )
        +sum( planth[gh, :node]==bus[n,:node] ? php[gh,t] : 0 for gh in NGH )
        + pshed[n,t] - pcc[n,t])
        == deltaT*sum( (exist(n,mapbranch[br,2],mapbranch[[br],:])==1) ? f[br,t] : 0 for br in NBR)  
           -deltaT*sum( (exist(mapbranch[br,1],n,mapbranch[[br],:])==1) ? f[br,t] : 0 for br in NBR)); 
#bolsa
if cons_efi=="Eficiencia"
    #@constraint(CPF,bolsa_gn[tt in ND], sum(plant[g,:Bolsa]=="B_Camisea" ? sum( tt==time_chart_df[t,:nday] ? cg[g,t] : 0 for t in NPER) : 0 for g in NGEN)<= bolsagn);
    @constraint(CPF,bolsa_gn_d[g in NGEN], sum(plant[g,:Bolsa]=="B_Camisea" ? cg[g,t] : 0 for t in NPER)<=cg_bolsa_d[g]);           
end
                       
#Cost           
@constraint(CPF,CostE==(deltaT*sum(plant[g,:cvnc]*pt[g,t] for g in NGEN,t in NPER)+sum(plant[g,:Co]*cg[g,t] for g in NGEN,t in NPER)
                    +sum(plant[g,:C_start]*y[g,t] for g in NGEN,t in NPER)+sum(plant[g,:C_down]*w[g,t] for g in NGEN,t in NPER)
                    +deltaT*sum(pshed[n,t]*price_shed for n in NBUS, t in NPER)
                    +deltaT*sum(planth[gh,:cost]*php[gh,t] for gh in NGH,t in NPER)));
#Objective
@objective(CPF, Min, CostE);
optimize!(CPF);
global lmpp=JuMP.dual.(pbalance);
global lfl1=JuMP.dual.(flowcap1);
global lfl2=JuMP.dual.(flowcap2);
global pshed_cx2=JuMP.value.(pshed);


global elapsed_time = Base.time_ns() - start_time
println("Elapsed time: ", elapsed_time / 1e9, " segundos")
include(string(folder_name,"/","print_acoplado.jl"))

#= # Escribir modelo en archivo de texto
open("modelo.txt", "w") do io
    print(io, PG)
end =#