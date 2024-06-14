       #Results of power generation
        pt_results=Array{Any}(nothing,nper+1,ngen+1);
        pt_results[1,1]=string("periods");
        for g in NGEN, t in NPER
                pt_results[1,g+1]=string(plant[g, :Unit]);
                pt_results[t+1,[1 g+1]]=[t pt_cx[g,t]];
        end
        pt_results_df = DataFrame(round.(pt_results[2:end,:], digits=digits_round), :auto);
        rename!( pt_results_df, Symbol.(pt_results[1,:]));
        CSV.write(string(path_output,"Thermic Generation(MW).csv"), pt_results_df);
        
        #Results of vol
        ngem=size(hiemb,1)
if ngem>=1
        vol_results=Array{Any}(nothing,nper+1,ngem+1);
        vol_results[1,1]=string("periods");
        for g in 1:ngem, t in NPER
                vol_results[1,g+1]=string( hiemb[g, :Unit]);
                vol_results[t+1,[1 g+1]]=[t vol_cx[NGHE[g],t]];
        end
        vol_results_df = DataFrame(round.(vol_results[2:end,:], digits=digits_round), :auto);
        rename!( vol_results_df, Symbol.(vol_results[1,:]));
        CSV.write(string(path_output,"Volumen(MMm3).csv"), vol_results_df);

else
        vol_results=Array{Any}(nothing,nper+1,ngh+1);
        vol_results[1,1]=string("periods");
        for g in 1:ngh, t in NPER
                vol_results[1,g+1]=string( planth[g, :Unit]);
                vol_results[t+1,[1 g+1]]=[t planth[g, :Vmax]];
        end
        vol_results_df = DataFrame(round.(vol_results[2:end,:], digits=digits_round), :auto);
        rename!( vol_results_df, Symbol.(vol_results[1,:]));
        CSV.write(string(path_output,"Volumen(MMm3).csv"), vol_results_df);        
end

if nrer>=1
        #Results of rer generation
        prer_results=Array{Any}(nothing,nper+1,nrer+1);
        prer_results[1,1]=string("periods");
        for g in NRER, t in NPER
                prer_results[1,g+1]=string( rer[g, :unit]);
                prer_results[t+1,[1 g+1]]=[t prer_cx[g,t]];
        end
        prer_results_df = DataFrame(round.(prer_results[2:end,:], digits=digits_round), :auto);
        rename!( prer_results_df, Symbol.(prer_results[1,:]));
        CSV.write(string(path_output,"RER generation(MW).csv"), prer_results_df);
end

if ngh>=1
        #Results of ghp generation
        pghp_results=Array{Any}(nothing,nper+1,ngh+1);
        pghp_results[1,1]=string("periods");
        for g in NGH, t in NPER
               pghp_results[1,g+1]=string(planth[g, :Unit]);
               pghp_results[t+1,[1 g+1]]=[t php_cx[g,t]];
        end
        pghp_results_df = DataFrame(round.(pghp_results[2:end,:], digits=digits_round), :auto);
        rename!( pghp_results_df, Symbol.(pghp_results[1,:]));
        CSV.write(string(path_output,"Hidralic Generation(MW).csv"), pghp_results_df);
end

        #Results of pressures
        pre_results=Array{Any}(nothing,nper+1,ngnode+1);
        pre_results[1,1]=string("periods");
        for gn in NGN, t in NPER
                pre_results[1,gn+1]=string(gnode[gn, :Node]);
                pre_results[t+1,[1 gn+1]]=[t pr_cx[gn,t]/fac_convert];
        end
        pre_results_df = DataFrame(round.(pre_results[2:end,:], digits=digits_round+1), :auto);
        rename!( pre_results_df, Symbol.(pre_results[1,:]));
        CSV.write(string(path_output,"Pressure (bar).csv"), pre_results_df);

        #Results of supply gas
        sup_results=Array{Any}(nothing,nper+1,nwell+1);
        sup_results[1,1]=string("periods");
        for s in NWE, t in NPER
                sup_results[1,s+1]=string(well[s, :Well]);
                sup_results[t+1,[1 s+1]]=[t su_cx[s,t]];
        end
        sup_results
        sup_results_df = DataFrame(round.(sup_results[2:end,:], digits=digits_round), :auto)
        rename!( sup_results_df, Symbol.(sup_results[1,:]))
        CSV.write(string(path_output,"Supply of gas(Mm3).csv"), sup_results_df)

        #Consum of power generation
        cg_results=Array{Any}(nothing,nper+1,ngen+1);
        cg_results[1,1]=string("periods");
        for g in NGEN, t in NPER
                cg_results[1,g+1]=string( plant[g, :Unit]);
                cg_results[t+1,[1 g+1]]=[t cg_cx[g,t]];
        end
        cg_results_df = DataFrame(round.(cg_results[2:end,:], digits=digits_round+1), :auto);
        rename!( cg_results_df, Symbol.(cg_results[1,:]));
        CSV.write(string(path_output,"Power consum of gas(Mm3).csv"), cg_results_df);

        #Results of start-up
        pt_start_up_results=Array{Any}(nothing,nper+1,ngen+1);
        pt_start_up_results[1,1]=string("periods");
        for g in NGEN, t in NPER
                pt_start_up_results[1,g+1]=string("gen ", plant[g, :Unit]);
                pt_start_up_results[t+1,[1 g+1]]=[t y_cx[g,t]];
        end
        pt_start_up_results_df = DataFrame(round.(pt_start_up_results[2:end,:], digits=digits_round-2), :auto);
        rename!( pt_start_up_results_df, Symbol.(pt_start_up_results[1,:]));
        CSV.write(string(path_output,"Start-up.csv"), pt_start_up_results_df);
        
        #Results of shut-down
        pt_shut_down_up_results=Array{Any}(nothing,nper+1,ngen+1);
        pt_shut_down_up_results[1,1]=string("periods");
        for g in NGEN, t in NPER
                pt_shut_down_up_results[1,g+1]=string("gen ", plant[g, :Unit]);
                pt_shut_down_up_results[t+1,[1 g+1]]=[t w_cx[g,t]];
        end
        pt_shut_down_up_results_df = DataFrame(round.(pt_shut_down_up_results[2:end,:], digits=digits_round-2), :auto);
        rename!( pt_shut_down_up_results_df, Symbol.(pt_shut_down_up_results[1,:]));
        CSV.write(string(path_output,"Shut_down.csv"), pt_shut_down_up_results_df);

        #Results of gas shed
        gshed_results=Array{Any}(nothing,nper+1,ngnode+1)
        gshed_results[1,1]=string("periods")
        for gn in NGN, t in NPER
                gshed_results[1,gn+1]=string("gnode ", gnode[gn, :Node])
                gshed_results[t+1,[1 gn+1]]=[t gshed_cx[gn,t]]
        end
        gshed_results_df = DataFrame(round.(gshed_results[2:end,:], digits=digits_round), :auto);
        rename!( gshed_results_df, Symbol.(gshed_results[1,:]));
        CSV.write(string(path_output,"Gas shed(MM3).csv"), gshed_results_df);

        #Results of power shed
        pshed_results=Array{Any}(nothing,nper+1,nbus+1);
        pshed_results[1,1]=string("periods");
        for b in NBUS, t in NPER
                pshed_results[1,b+1]=string("bus ", bus[b, :node]);
                pshed_results[t+1,[1 b+1]]=[t pshed_cx[b,t]];
        end
        pshed_results_df = DataFrame(round.(pshed_results[2:end,:], digits=digits_round), :auto);
        rename!( pshed_results_df, Symbol.(pshed_results[1,:]));
        CSV.write(string(path_output,"Power shed(MW).csv"), pshed_results_df);

        #Results of power flow
        flowdc_results=Array{Any}(nothing,nper+1,nbranch+1);
        flowdc_results[1,1]=string("periods");
        for br in NBR,t in NPER
            flowdc_results[1,br+1]=string(branch[br,:branch]," ",mapbranch[br,1],"->",mapbranch[br,2]);
            flowdc_results[t+1,[1 br+1]]=[t f_cx[br,t]];
        end
        flowdc_results_df = DataFrame(round.(flowdc_results[2:end,:], digits=digits_round), :auto);
        rename!(flowdc_results_df, Symbol.(flowdc_results[1,:]));
        CSV.write(string(path_output,"DC power flow branches(MW).csv"), flowdc_results_df);
        
        #Results of gas flow
        flowg_results=Array{Any}(nothing,nper+1,npipe+1);
        flowg_results[1,1]=string("periods");
        for pp in NPP, t in NPER
            flowg_results[1,pp+1]=string(pipe[pp,:pipe]," ",mappipe[pp,1],"->",mappipe[pp,2]);
            flowg_results[t+1,[1 pp+1]]=[t fg_cx[pp,t]];
        end
        flowg_results_df = DataFrame(round.(flowg_results[2:end,:], digits=digits_round+1), :auto);
        rename!(flowg_results_df, Symbol.(flowg_results[1,:]));
        CSV.write(string(path_output,"Gas flow pipe(Mm3).csv"), flowg_results_df);

        #Results of line pack
        lpg_results=Array{Any}(nothing,nper+1,npipe+1);
        lpg_results[1,1]=string("periods");
        for pp in NPP, t in NPER
            lpg_results[1,pp+1]=string(pipe[pp, :pipe]);
            lpg_results[t+1,[1 pp+1]]=[t lip_cx[pp,t]];
        end
        lpg_results_df = DataFrame(round.(lpg_results[2:end,:], digits=digits_round+1), :auto);
        rename!(lpg_results_df, Symbol.(lpg_results[1,:]));
        CSV.write(string(path_output,"Line pack (Mm3).csv"), lpg_results_df);
        

        fig_results=Array{Any}(nothing,nper+1,npipe+1);
        fig_results[1,1]=string("periods");
        for pp in NPP, t in NPER
            fig_results[1,pp+1]=string("pipe ", pipe[pp, :pipe]);
            fig_results[t+1,[1 pp+1]]=[t fig_cx[pp,t]];
        end
        fig_results_df = DataFrame(round.(fig_results[2:end,:],digits=digits_round+1), :auto)
        rename!(fig_results_df, Symbol.(fig_results[1,:]));
        CSV.write(string(path_output,"Gas flow in (Mm3).csv"), fig_results_df);

if ncomp>=1
        #Flujo de compresor
        fco_results=Array{Any}(nothing,nper+1,ncomp+1);
        fco_results[1,1]=string("periods");
        for com in NCO, t in NPER
            fco_results[1,com+1]=string("pipe ", comp[com, :comp]);
            fco_results[t+1,[1 com+1]]=[t fco_cx[com,t]];
        end
        fco_results_df = DataFrame(round.(fco_results[2:end,:],digits=digits_round+1), :auto)
        rename!(fco_results_df, Symbol.(fco_results[1,:]))
        CSV.write(string(path_output,"Flujo de compresor (Mm3).csv"), fco_results_df);

        #Consumo de compresor
        cco_results=Array{Any}(nothing,nper+1,ncomp+1);
        cco_results[1,1]=string("periods");
        for com in NCO, t in NPER
            cco_results[1,com+1]=string("pipe ", comp[com, :comp]);
            cco_results[t+1,[1 com+1]]=[t cco_cx[com,t]];
        end
        cco_results_df = DataFrame(round.(cco_results[2:end,:],digits=digits_round+1), :auto)
        rename!(cco_results_df, Symbol.(cco_results[1,:]))
        CSV.write(string(path_output,"Consumo de compresor (Mm3).csv"), cco_results_df);
    end

if nreg>=1
        #Flujo de regulador
        fre_results=Array{Any}(nothing,nper+1,nreg+1);
        fre_results[1,1]=string("periods");
        for re in NRE, t in NPER
            fre_results[1,re+1]=string("pipe ", reg[re, :rep]);
            fre_results[t+1,[1 re+1]]=[t fre_cx[re,t]];
        end
        fre_results_df = DataFrame(round.(fre_results[2:end,:],digits=digits_round+1), :auto)
        rename!(fre_results_df, Symbol.(fre_results[1,:]))
        CSV.write(string(path_output,"Flujo de Regulador (Mm3).csv"), fre_results_df);

        #Flujo de regulador
        cre_results=Array{Any}(nothing,nper+1,nreg+1);
        cre_results[1,1]=string("periods");
        for re in NRE, t in NPER
            cre_results[1,re+1]=string("pipe ", reg[re, :rep]);
            cre_results[t+1,[1 re+1]]=[t cre_cx[re,t]];
        end
        cre_results_df = DataFrame(round.(cre_results[2:end,:],digits=digits_round+1), :auto)
        rename!(cre_results_df, Symbol.(cre_results[1,:]))
        CSV.write(string(path_output,"Consum de Regulador (Mm3).csv"), cre_results_df);
end

        #Results of flujo de salida de gas
        fog_results=Array{Any}(nothing,nper+1,npipe+1);
        fog_results[1,1]=string("periods");
        for pp in NPP, t in NPER
            fog_results[1,pp+1]=string("pipe ", pipe[pp, :pipe]);
            fog_results[t+1,[1 pp+1]]=[t fog_cx[pp,t]];
        end
        fog_results_df = DataFrame(round.(fog_results[2:end,:],digits=digits_round+1), :auto);
        rename!(fog_results_df, Symbol.(fog_results[1,:]));
        CSV.write(string(path_output,"Gas flow out (Mm3).csv"), fog_results_df);

       #Results of LMP (Active power)
       lmpp_results=Array{Any}(nothing,nper+1,nbus+1);
       lmpp_results[1,1]=string("periods");
       for n in NBUS, t in NPER
                lmpp_results[1,n+1]=string("bus ", bus[n, :node]);
                lmpp_results[t+1,[1 n+1]]=[t lmpp[n,t]];
        end
        lmpp_results_df = DataFrame(round.(lmpp_results[2:end,:], digits=digits_round), :auto);
        rename!(lmpp_results_df, Symbol.(lmpp_results[1,:]));
        CSV.write(string(path_output,"LMP(Dolar per MWh).csv"), lmpp_results_df);

        #Results of flow cap 1
        lmp1_results=Array{Any}(nothing,nper+1,nbranch+1);
        lmp1_results[1,1]=string("periods");
        for n in NBR, t in NPER
                lmp1_results[1,n+1]=string("id",branch[n,:branch]," ",mapbranch[n,1],"->",mapbranch[n,2]);
                lmp1_results[t+1,[1 n+1]]=[t lfl1[n,t]];
        end
        lmp1_results_df = DataFrame(round.(lmp1_results[2:end,:], digits=digits_round), :auto);
        rename!(lmp1_results_df, Symbol.(lmp1_results[1,:]));
        CSV.write(string(path_output,"LMP Flowcap_1.csv"), lmp1_results_df);
   
        #Results of flow cap 2
        lmp2_results=Array{Any}(nothing,nper+1,nbranch+1);
        lmp2_results[1,1]=string("periods");
        for n in NBR, t in NPER
                lmp2_results[1,n+1]=string("id",branch[n,:branch]," ",mapbranch[n,1],"->",mapbranch[n,2]);
                lmp2_results[t+1,[1 n+1]]=[t lfl2[n,t]];
        end
        lmp2_results_df = DataFrame(round.(lmp2_results[2:end,:], digits=digits_round), :auto);
        rename!(lmp2_results_df, Symbol.(lmp2_results[1,:]));
        CSV.write(string(path_output,"LMP Flowcap_2.csv"), lmp2_results_df);

        viol=Array{Any}(nothing,npipe,nper);
        for pp in NPP, t in NPER
            viol[pp,t]=[abs((pr_cx[mappipe[pp,1],t])^2 - (pipe[pp,:sqrt_es]*pr_cx[mappipe[pp,2],t])^2-(fg_cx[pp,t])^2/(pipe[pp,:Cmn]/stageG)^2)/((pr_cx[mappipe[pp,1],t])^2 )];
        end
        Violmax=maximum(viol,dims=2)
        ViolmaxT=maximum(Violmax,dims=1)[1]      
        Violmean1=sum(viol)/nper/npipe
        
        #Violation of pipeline constraint
        del_results=Array{Any}(nothing,nper+1,npipe+1);
        del_results[1,1]=string("periods");
        for pp in NPP, t in NPER
            del_results[1,pp+1]=string("pipe ", pipe[pp, :pipe]);
            del_results[t+1,[1 pp+1]]=[t viol[pp,t]];
        end
        del_results_df = DataFrame(round.(del_results[2:end,:],digits=digits_round+1), :auto);
        rename!(del_results_df, Symbol.(del_results[1,:]));
        CSV.write(string(path_output,"Violation of pipeline constraint.csv"), del_results_df);

       #Results of LMP (GN)
       lmpp_results=Array{Any}(nothing,nper+1,ngnode+1);
       lmpp_results[1,1]=string("periods");
       for n in NGN, t in NPER
                lmpp_results[1,n+1]=string( gnode[n, :Node]);
                lmpp_results[t+1,[1 n+1]]=[t lmpp_gn[n,t]];
        end
        lmpp_results_df = DataFrame(round.(lmpp_results[2:end,:], digits=digits_round+5), :auto);
        rename!(lmpp_results_df, Symbol.(lmpp_results[1,:]));
        CSV.write(string(path_output,"LMP(Dolar per M3).csv"), lmpp_results_df);

        #Results of flow cap 1
        lmp1_results=Array{Any}(nothing,nper+1,ngnode+1);
        lmp1_results[1,1]=string("periods");
        for n in NGN, t in NPER
                lmp1_results[1,n+1]=string("ngnode ", gnode[n, :Node]);
                lmp1_results[t+1,[1 n+1]]=[t lp_min_gn[n,t]];
        end
        lmp1_results_df = DataFrame(round.(lmp1_results[2:end,:], digits=digits_round), :auto);
        rename!(lmp1_results_df, Symbol.(lmp1_results[1,:]));
        CSV.write(string(path_output,"LMP Pressure min.csv"), lmp1_results_df);
   
        #Results of flow cap 2
        lmp2_results=Array{Any}(nothing,nper+1,ngnode+1);
        lmp2_results[1,1]=string("periods");
        for n in NGN, t in NPER
                lmp2_results[1,n+1]=string("ngnode ", gnode[n, :Node]);
                lmp2_results[t+1,[1 n+1]]=[t lp_max_gn[n,t]];
        end
        lmp2_results_df = DataFrame(round.(lmp2_results[2:end,:], digits=digits_round), :auto);
        rename!(lmp2_results_df, Symbol.(lmp2_results[1,:]));
        CSV.write(string(path_output,"LMP Pressure max.csv"), lmp2_results_df);

        #Results of Objective Function
        OF_results=Array{Any}(nothing,10,2);
        OF_results[1,1]="Parameter";
        OF_results[1,2]="Value";
        OF_results[2,1]="Power shedding cost MM (USD)";
        OF_results[2,2]=round.(price_shed*deltaT*sum(pshed_cx[n,t] for n in NBUS, t in NPER)/1000000,digits=3);
        OF_results[3,1]="Thermal variable cost MM (USD)";
        OF_results[3,2]=round.(deltaT*sum(plant[g,:cvnc]*pt_cx[g,t] for g in NGEN,t in NPER)/1000000+deltaT*sum(plant[g,:Co]*cg_cx[g,t] for g in NGEN,t in NPER)/1000000,digits=3);
        OF_results[4,1]="Start-up and-Shut down MM (USD)";
        OF_results[4,2]=round.((sum(plant[g,:C_start]*y_cx[g,t] for g in NGEN,t in NPER)+sum(plant[g,:C_down]*w_cx[g,t] for g in NGEN,t in NPER))/1000000,digits=3);
        OF_results[5,1]="Gas shedding cost MM (USD)";
        OF_results[5,2]=round.(g_price_shed*sum(gshed_cx[gn,t] for gn in NGN, t in NPER)/1000000,digits=3);
        OF_results[6,1]="Well Production MM (USD)";
        OF_results[6,2]=round.(sum(well[w,:C_sup_D]*su_cx[w,t] for w in NWE,t in NPER)/1000000,digits=3);
        OF_results[7,1]="Hidraulic Production MM (USD)";
        OF_results[7,2]=if ngh+ngem>=1 round.(sum(planth[gh,:cost]*php_cx[gh,t] for gh in NGH,t in NPER)/1000000,digits=3) else 0 end;
        OF_results[8,1]="Total cost MM (USD)";
        OF_results[8,2]=OF_results[2,2]+OF_results[3,2] +OF_results[4,2] +OF_results[5,2] +OF_results[6,2] +OF_results[7,2] 
        OF_results[9,1]="Time seg";
        OF_results[9,2]=elapsed_time/1e9; 
        OF_results[10,1]="% viol conx gas flow max";
        OF_results[10,2]=ViolmaxT[1]*1e2; 
        OF_results_df = DataFrame(OF_results[2:end,:], :auto)
        rename!(OF_results_df, Symbol.(OF_results[1,:]));
        CSV.write(string(path_output,"Objective function.csv"), OF_results_df);

#Results of Cost per time
subcost_names=["Hidraulic Power" "Thermal power" "Start-up and-Shut down" "Power shedding" "Well Production" "Gas shedding"]
        subcost = 6

        cost_t_0=Array{Any}(nothing,subcost,nper);
        for t in NPER
            cost_t_0[1,t]=[if ngh + ngem >= 1 round(sum(planth[gh, :cost] * php_cx[gh, t] for gh in NGH) / 1000000, digits = 6) else 0 end];
            cost_t_0[2,t]=[round(deltaT * sum(plant[g, :cvnc] * pt_cx[g, t] for g in NGEN) / 1000000 + deltaT * sum(plant[g, :Co] * cg_cx[g, t] for g in NGEN) / 1000000, digits = 6)];
            cost_t_0[3,t]=[round(sum(plant[g, :C_start] * y_cx[g, t] for g in NGEN)/1000000 + sum(plant[g, :C_down] * w_cx[g, t] for g in NGEN) / 1000000, digits = 6)];
            cost_t_0[4,t]=[round(price_shed * deltaT * sum(pshed_cx[n, t] for n in NBUS) / 1000000, digits = 6)];
            cost_t_0[5,t]=[round(sum(well[w, :C_sup_D] * su_cx[w, t] for w in NWE) / 1000000, digits = 6)];
            cost_t_0[6,t]=[round(g_price_shed * sum(gshed_cx[gn, t] for gn in NGN) / 1000000, digits = 6)];
        end

        Cost_h = Array{Any}(nothing, nper + 1, subcost + 1)
        Cost_h[1, 1] = "periods"
        for sc in 1:subcost, t in NPER
                Cost_h[1, sc + 1] = subcost_names[sc]
                Cost_h[t + 1, [1 sc+1]] = [t cost_t_0[sc,t]]
        end     
        Cost_h_df = DataFrame(Cost_h[2:end, :], :auto);
        rename!(Cost_h_df, Symbol.(Cost_h[1, :]));
        CSV.write(string(path_output, "Cost_per_t.csv"), Cost_h_df);

if folder_end=="Acoplado_eMISOCP"
        #Results of iters
        ite_results=Array{Any}(nothing,iter+1,4);
        ite_results[1,1]=string("iterations");
        ite_results[1,2]=string("Cost MMUSD");
        ite_results[1,3]=string("% Error Mean");
        ite_results[1,4]=string("% Error Max");
        for j in 1:iter
            ite_results[j+1,1]=j;
            ite_results[j+1,2]=costs[j]/10^6;
            ite_results[j+1,3]=error_mean[j]*100;
            ite_results[j+1,4]=error_max[j][1]*100;
        end
        ite_results_df = DataFrame(round.(ite_results[2:end,:], digits=digits_round+1), :auto);
        rename!( ite_results_df, Symbol.(ite_results[1,:]));
        CSV.write(string(path_output,"Iterations.csv"), ite_results_df);
end

       #Apoyo
       ap_results=Array{Any}(nothing,nper+1,4);
       lp_total=Array{Any}(nothing,nper);
       ap_results[1,1]=string("periods");
       lp_total=vec(sum(lip_cx.data, dims=1))
       for t in NPER
                ap_results[1,3]="ConsumSP";
                ap_results[1,2]="Delta LP";
                ap_results[1,4]="Delta Volumen";
                ap_results[t+1,1]=t;
                ap_results[t+1,3]=sum(plant[g,:Bolsa]=="B_Camisea" ? cg_cx[g,t] : 0 for g in NGEN);
                ap_results[t+1,2]=lp_total[t] - ( t==1 ? sum(pipe[:,:L_0], dims=1)[1] : lp_total[t-1]);
                if ngem>1
                        vol_total=vec(sum(vol_cx.data, dims=1))
                        ap_results[t+1,4]=vol_total[t] - ( t==1 ? sum(planth[:,:Vi], dims=1)[1] : vol_total[t-1])
                else
                        ap_results[t+1,4]=0
                end
       end
       ap_results_df = DataFrame(round.(ap_results[2:end,:], digits=digits_round), :auto);
       rename!( ap_results_df, Symbol.(ap_results[1,:]));
       CSV.write(string(path_output,"Apoyo.csv"), ap_results_df);

               #Results of LMP (Active power)
               lmpp_results=Array{Any}(nothing,nper+1,nbus+1);
               lmpp_results[1,1]=string("periods");
               for n in NBUS, t in NPER
                        lmpp_results[1,n+1]=string("bus ", bus[n, :node]);
                        lmpp_results[t+1,[1 n+1]]=[t lmpp_a[n,t]];
                end
                lmpp_results_df = DataFrame(round.(lmpp_results[2:end,:], digits=digits_round), :auto);
                rename!(lmpp_results_df, Symbol.(lmpp_results[1,:]));
                CSV.write(string(path_output,"LMP Acoplado(Dolar per MWh).csv"), lmpp_results_df);
                
        #Results of gas flow
        lp_flowg_results=Array{Any}(nothing,nper+1,npipe+1);
        lp_flowg_results[1,1]=string("periods");
        for pp in NPP, t in NPER
            lp_flowg_results[1,pp+1]=string(pipe[pp,:pipe]," ",mappipe[pp,1],"->",mappipe[pp,2]);
            lp_flowg_results[t+1,[1 pp+1]]=[t lp_flow_gn[pp,t]];
        end
        lp_flowg_results_df = DataFrame(round.(lp_flowg_results[2:end,:], digits=digits_round+1), :auto);
        rename!(lp_flowg_results_df, Symbol.(lp_flowg_results[1,:]));
        CSV.write(string(path_output,"LMP Gas flow pipe.csv"), lp_flowg_results_df);