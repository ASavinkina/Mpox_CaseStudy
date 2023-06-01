#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(adaptivetau)
library(ggplot2)
library(tidyverse)
library(shinydashboard)
library(scales)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    output$I_plot <- renderPlot({
        
        ############# Multiple plot function ###################
        ## http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/ ##
        
        multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
            library(grid)
            
            # Make a list from the ... arguments and plotlist
            plots <- c(list(...), plotlist)
            
            numPlots = length(plots)
            
            # If layout is NULL, then use 'cols' to determine layout
            if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
            }
            
            if (numPlots==1) {
                print(plots[[1]])
                
            } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
                
                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                    # Get the i,j matrix positions of the regions that contain this subplot
                    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                    
                    print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                    layout.pos.col = matchidx$col))
                }
            } 
        }
    
    
    # Initial conditions
    
    Yale_undegrad_pop = input$Yale_undegrad_pop #input$Yale_undegrad_pop # total population
    initial_inf= input$initial_inf #input$initial_inf # initial infections
    Vaccine_Efficacy = input$vaxefficacy
    Vaccine_Efficacy_PEP = Vaccine_Efficacy * 0.5
    exogshock =1 #size of exogenous shock. 0 for no exogenous shocks, >1 for superspreader events
    exograte = ifelse(input$exograte==0, 0, input$exograte/100)  # rate of exogenous shocks
    recoveryrate = 1/21
    diagrate = ifelse(input$diagrate>0.95, 1, (input$diagrate*recoveryrate)/(1-input$diagrate)) # daily rate of diagnosis of infectious cases
    isorate = 1/input$isoduration #duration of isolation
    studentcontacts = input$studentcontacts # students vaccinated per diagnosed case
    R0_h= input$R0_h
    R0_l= input$R0_l
    HR_prop = input$HR_prop
    LR_prop = 1-input$HR_prop
    Prop_Vaccinated = input$propVax
    Cost_vax = input$costvax
    Cost_detect = input$costdetect
      

    init.values = c(
      S_h = Yale_undegrad_pop*HR_prop-initial_inf-Yale_undegrad_pop*HR_prop*Vaccine_Efficacy*Prop_Vaccinated, P_h = 0, I_h = initial_inf, Dx0_h=0, Dx_h=0, 
      Vs_h=Yale_undegrad_pop*HR_prop*Vaccine_Efficacy*Prop_Vaccinated, Vi_h=0, R_h=0,
      S_l = Yale_undegrad_pop*LR_prop, P_l = 0, I_l = 0, Dx0_l=0, Dx_l=0, 
      Vs_l=0, Vi_l=0, R_l=0
    )
    
    # Specify all transitions
    transitions = list(
      # high risk population
      c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, endogenous infection
      c(S_h = -1, P_h = +1), # movement from susceptible to presymptomatic, exogenous infection 
      c(S_h= -1, Vs_h= +1), #movement from susceptible to vaccinated susceptible
      c(P_h= -1, Vi_h= +1), #movement from presymptomatuc to vaccinated infected
      c(P_h = -1, I_h = +1), #movement from presymptomatic to infected 
      c(I_h = -1, Dx0_h = +1),  #movement from infected to newly diagnosed
      c(Dx0_h = -1, Dx_h = +1), #movement from newly diagnosed to diagnosed
      c(I_h = -1, R_h = +1), #movement from infected to recovered
      c(Dx_h= -1, R_h= +1), #movement from diagnosed to recovered
      
      #low risk population
      c(S_l = -1, P_l = +1), # movement from susceptible to presymptomatic, endogenous infection
      c(S_l= -1, Vs_l= +1), #movement from susceptible to vaccinated susceptible
      c(P_l= -1, Vi_l= +1), #movement from presymptomatuc to vaccinated infected
      c(P_l = -1, I_l = +1), #movement from presymptomatic to infected 
      c(I_l = -1, Dx0_l = +1),  #movement from infected to newly diagnosed
      c(Dx0_l = -1, Dx_l = +1), #movement from newly diagnosed to diagnosed
      c(I_l = -1, R_l = +1), #movement from infected to recovered
      c(Dx_l= -1, R_l= +1) #movement from diagnosed to recovered
        
    )
    
    # Rates for all transitions
    # (In same order as "transitions")
    RateF <- function(x, pars, times) {
      return(c(
        pars$beta_l*x["S_h"]*x["I_l"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"])+
        pars$beta_l*x["S_h"]*x["I_h"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"])+
        pars$beta_h*x["S_h"]*x["I_h"]/(x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"]),                                                                                
        pars$theta, # movement from susceptible to presymptomatic, exogenous infection
        ifelse(x["P_h"]>1,pars$iota*x["Dx0_h"]*(1-pars$attackrate)*pars$Vaccine_Efficacy,pars$iota*x["Dx0_h"]*pars$Vaccine_Efficacy),
        ifelse(x["P_h"]>1, pars$iota*x["Dx0_h"]*pars$attackrate*pars$Vaccine_Efficacy_PEP, 0), #movement from susceptible to vaccinated infected (based on number newly diagnosed)
        pars$gamma*x["P_h"], #movement from presymptomatic to infected (duration of incubation)
        pars$delta*x["I_h"], #movement from infected to diagnosed
        pars$tau*x["Dx0_h"], #movement from newly diagnosed to diagnosed
        pars$rho*x["I_h"], #movement from infected to recovered
        pars$omicron*x["Dx_h"], #movement from diagnosed to recovered 
        
        pars$beta_l*x["S_l"]*x["I_l"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"])+
          pars$beta_l*x["S_l"]*x["I_h"]/(x["S_l"] + x["P_l"] + x["I_l"] + x["Vi_l"]+  x["Vs_l"] + x["R_l"]+ x["S_h"] + x["P_h"] + x["I_h"] + x["Vi_h"]+  x["Vs_h"]+ x["R_h"]),
        ifelse(x["P_l"]>1,pars$iota*x["Dx0_l"]*(1-pars$attackrate)*pars$Vaccine_Efficacy,pars$iota*x["Dx0_l"]*pars$Vaccine_Efficacy),
        ifelse(x["P_l"]>1, pars$iota*x["Dx0_l"]*pars$attackrate*pars$Vaccine_Efficacy_PEP, 0), #movement from susceptible to vaccinated infected (based on number newly diagnosed)
        pars$gamma*x["P_l"], #movement from presymptomatic to infected (duration of incubation)
        pars$delta*x["I_l"], #movement from infected to diagnosed
        pars$tau*x["Dx0_l"], #movement from newly diagnosed to diagnosed
        pars$rho*x["I_l"], #movement from infected to recovered
        pars$omicron*x["Dx_l"] #movement from diagnosed to recovered 
      ))
    }
    
    # Setting parameters
    pars = list(
      beta_h= R0_h * (1/21),  # beta
      beta_l = R0_l * (1/21),
      gamma= 1/7.6, #incubation period
      delta= diagrate,  #diagnosis rate
      rho= recoveryrate, #recovery rate
      tau = 1, #rate of moving from newly diagnosed to diagnosed
      iota= studentcontacts, # students vaccinated per diagnosed case
      mu = 1/7.6, #incubation period for those in vaccination
      #omega = quarrate, #length of vaccination for susceptible
      attackrate=0.2, 
      theta= exograte, #rate of exogenous shocks
      omicron = isorate,
      Vaccine_Efficacy=Vaccine_Efficacy,
      Vaccine_Efficacy_PEP = Vaccine_Efficacy_PEP
    )
    
    
    # Running stochastic model
    
    
    # results = as.data.frame(ssa.adaptivetau(init.values, 
    #                                         transitions, 
    #                                         RateF, 
    #                                         pars, 
    #                                         tf=100))
    
    #Running stochastic model 1,000 times
    
    #  Create dataset
    
    
    results_all_presymptomatic <- data.frame(matrix(0, nrow=0, ncol=4))
    colnames(results_all_presymptomatic) <- c("time", "P_h","P_l", "run")
    results_all_infected <- data.frame(matrix(0, nrow=0, ncol=4))
    colnames(results_all_infected) <- c("time", "I_h","I_l", "run")
    results_all_diagnosed <- data.frame(matrix(0, nrow=0, ncol=4))
    colnames(results_all_diagnosed) <- c("time", "Dx_h","Dx_l", "run")
    results_all_recovered <- data.frame(matrix(0, nrow=0, ncol=4))
    colnames(results_all_recovered) <- c("time", "R_h","R_l", "run")
    results_all_newlydiagnosed <- data.frame(matrix(0, nrow=0, ncol=4))
    colnames(results_all_newlydiagnosed) <- c("time", "Dx0_h","Dx0_l", "run")
    results_all_vaccinated<- data.frame(matrix(0, nrow=0, ncol=4))
    colnames(results_all_vaccinated) <- c("time", "Vx_h","Vx_l", "run")
    
    runs=100
    
    set.seed(100)
    
    
    for (i in 1:runs) {
        
        
        results = as.data.frame(ssa.adaptivetau(init.values, 
                                                transitions, 
                                                RateF, 
                                                pars, 
                                                tf=100))
        
        results_presymptomatic_i <- results[,c(1,3,11)]
        results_presymptomatic_i$run <- paste0(i)
        
        results_all_presymptomatic <- rbind(results_all_presymptomatic, results_presymptomatic_i)
        
        results_infected_i <- results[,c(1,4,12)]
        results_infected_i$run <- paste0(i)
        
        results_all_infected <- rbind(results_all_infected, results_infected_i)
        
        results_newlydiagnosed_i <- results[,c(1,5,13)]
        results_newlydiagnosed_i$run <- paste0(i)
        
        results_all_newlydiagnosed <- rbind(results_all_newlydiagnosed, results_newlydiagnosed_i)
        
        results_diagnosed_i <- results[,c(1,6,14)]
        results_diagnosed_i$run <- paste0(i)
        
        results_all_diagnosed <- rbind(results_all_diagnosed, results_diagnosed_i)
        
        results_recovered_i <- results[,c(1,9,17)] 
        results_recovered_i$run <- paste0(i)
        
        results_all_recovered <- rbind(results_all_recovered, results_recovered_i)
        
        results_vaccinated_i <- results[,c(1,7,15)]
        results_vaccinated_i[,2] <- results_vaccinated_i[,2]/Vaccine_Efficacy + results[,8]/Vaccine_Efficacy_PEP 
        results_vaccinated_i[,3] <- results_vaccinated_i[,3]/Vaccine_Efficacy + results[,16]/Vaccine_Efficacy_PEP 
        results_vaccinated_i$run <- paste0(i)
        
        results_all_vaccinated <- rbind(results_all_vaccinated, results_vaccinated_i)
        
    }
    

    # Total infections
    
    total_infections <- results_all_recovered
    
    total_infections$allinfections_hr <- total_infections$R_h + results_all_diagnosed$Dx_h + 
      results_all_newlydiagnosed$Dx0_h +
      results_all_infected$I_h +
      results_all_presymptomatic$P_h 
    
    total_infections$allinfections_lr <- total_infections$R_l + 
      results_all_diagnosed$Dx_l + 
      results_all_newlydiagnosed$Dx0_l +
      results_all_infected$I_l+ 
      results_all_presymptomatic$P_l
    
    #   Calculate total number of infectious students 
    
    total_infections2 <- total_infections[which(total_infections$time==100),]
    
    total_infections2$total <- total_infections2$allinfections_hr + total_infections2$allinfections_lr
    
    output$maxinfectionsplot <- renderPlot(hist(total_infections2$total,
                                                main="",xlab="Cumulative number of cases at 100 days",
                                                ylab="Percent likelihood", breaks=30),
                                           width=300, height=300)
    
    
    #   Calculate number of expected index cases
    
    expected_infections <- input$exograte+input$initial_inf
    
    #   Calculate proportion of times there are more cases than the expected number
    
    additional_cases_hr <- length(total_infections2$allinfections_hr[total_infections2$allinfections_hr>expected_infections])
    additional_cases_lr <- length(total_infections2$allinfections_lr[total_infections2$allinfections_lr>0])
    
    #   Calculate mean number of additional cases given secondary case
    
    additional_cases_mean_hr <- mean(total_infections2$allinfections_hr[total_infections2$allinfections_hr>expected_infections])
    additional_cases_mean_lr <- mean(total_infections2$allinfections_lr[total_infections2$allinfections_lr>expected_infections])
    
    cases_mean <- mean((total_infections2$allinfections_lr + total_infections$allinfections_hr))
    
    # Calculate min and max of total number of infections (assuming secondary infections)
    
    min_infections_hr <- min(total_infections2$allinfections_hr[total_infections2$allinfections_hr>expected_infections])
    max_infections_hr <- max(total_infections2$allinfections_hr[total_infections2$allinfections_hr>expected_infections])
    
    min_infections_lr <- min(total_infections2$allinfections_lr[total_infections2$allinfections_lr>0])
    max_infections_lr <- max(total_infections2$allinfections_lr[total_infections2$allinfections_lr>0])
    

    
    
    # Average vaccinate beds
    
    results_all_vaccinated2 <- results_all_vaccinated
    
    results_all_vaccinated2$V <- results_all_vaccinated2$Vs_h + results_all_vaccinated2$Vs_l
    
    results_all_vaccinated2$time2 <- round(results_all_vaccinated2$time)
    
    results_all_vaccinated3 <- results_all_vaccinated2 %>%
      group_by(run,time2) %>%
      summarise_at(vars(V), list(maxQ = max))
    
  
    # Likelihood exceeding vaccinate capacity
    
    results_all_vaccinated_likelihood <- results_all_vaccinated3 %>%
      group_by(run) %>%
      summarise_at(vars(maxQ), list(maxQQ = max))
    
    results_all_vaccinated_cost <- results_all_vaccinated_likelihood
    
    results_all_vaccinated_cost$cost <- results_all_vaccinated_likelihood$maxQQ* Cost_vax
    
    
    output$maxvaxplot <- renderPlot(hist(results_all_vaccinated_likelihood$maxQQ,
                                         main="",xlab="Max number of students vaccinated",
                                         ylab="Percent likelihood"), height=300, width=300)
    
             
    
    #likelihood_vax_past_cap <- percent(length(which(results_all_vaccinated_likelihood$maxQQ>vaccinate_capacity_count))/length(results_all_vaccinated_likelihood$maxQQ))
    
    results_all_vaccinated4 <- results_all_vaccinated3 %>%
      group_by(time2) %>%
      summarise_at(vars(maxQ), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                    median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                    max=max))

    avg_vaccinate_plot <- ggplot(data=results_all_vaccinated3, aes(x=time2, y=maxQ, group=run)) + geom_line() +
      theme_classic() + theme(legend.position = "none") + 
      #geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + 
      xlab("Days") + ylab("Average number of vaccinated students") +
      ggtitle("Number of vaccinated \nstudents by day")
    
    # Average infections students
    
    results_all_infected2 <- results_all_infected
    
    results_all_infected2$I <- results_all_infected$I_h + results_all_infected$I_l
    
    results_all_infected2$time2 <- round(results_all_infected2$time)
    
    results_all_infected3 <- results_all_infected2 %>%
      group_by(run,time2) %>%
      summarise_at(vars(I), list(maxI = max))
    
    results_all_infected4 <- results_all_infected3 %>%
      group_by(time2) %>%
      summarise_at(vars(maxI), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                    median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                    max=max))

    avg_infectious_plot <- ggplot(data=results_all_infected3, aes(x=time2, y=maxI, group=run)) + geom_line(alpha=0.3) +
      theme_classic() + theme(legend.position = "none") + 
      #geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + 
      xlab("Days") + ylab("Average number of infectious students") +
      ggtitle("Number of infectious \nstudents by day")
    
    
    # Average isolated students 
    
    #isolation_capacity_count <- input$isocapacity
    
    results_all_diagnosed2 <- cbind(results_all_diagnosed, results_all_newlydiagnosed$Dx0_h, results_all_newlydiagnosed$Dx0_l)
    
    results_all_diagnosed2$total <- results_all_diagnosed2$Dx_h + results_all_diagnosed2$Dx_l +
      results_all_diagnosed2$`results_all_newlydiagnosed$Dx0_h`+ 
      results_all_diagnosed2$`results_all_newlydiagnosed$Dx0_l`
    
    results_all_diagnosed2$time2 <- round(results_all_diagnosed2$time)
    
    results_all_diagnosed3 <- results_all_diagnosed2 %>%
      group_by(run,time2) %>%
      summarise_at(vars(total), list(maxD = max))
    
    results_all_diagnosed_likelihood <- results_all_diagnosed3 %>%
      group_by(run) %>%
      summarise_at(vars(maxD), list(maxDD = max))
    
  # Calculate cost of isolation
  # We are going to simplify this for the case study: isolation cost in model will be calculated as
  # # detected and isolated (ie number newly diagnosed) * cost of isolation per day * days of isolation
  # # This isn't quite accurate for the model since isolation time is a distribution, but it is close enough
    
    # Create a new dataset to calculate cost
    results_all_diagnosed_cost <- results_all_newlydiagnosed
    
    # Add together newly diagnosed high and low risk subgroups
    results_all_diagnosed_cost$newly_diagnosed <- results_all_diagnosed_cost$Dx0_h + results_all_diagnosed_cost$Dx0_l
    
    # Collapse time into days
    results_all_diagnosed_cost$time2 <- round(results_all_diagnosed_cost$time)
    
    # Keep maximum number isolated by day (some days have multiple events occur- multiple observations.
    # We want a daily count of number isolated. This won't be perfect but it will be close.)
    results_all_diagnosed_cost2 <- results_all_diagnosed_cost %>%
      group_by(run,time2) %>%
      summarise_at(vars(newly_diagnosed), list(newly_diag = max))
    
    # Calculate sum of people isolated for each run over entire time
    results_all_diagnosed_cost3 <- results_all_diagnosed_cost2 %>%
      group_by(run) %>%
      summarise(sum(newly_diag))
    
    # Calculate cost by multiplying number isolated * cost of isolation * time of isolation
    results_all_diagnosed_cost3$cost <- results_all_diagnosed_cost3$`sum(newly_diag)` * Cost_detect * input$isoduration
    
    

    output$maxisoplot <- renderPlot(hist(results_all_diagnosed_likelihood$maxDD,
                                         main="",xlab="Max number of students in isolation",
                                         ylab="Percent likelihood"), height=300, width=300)


    
    #likelihood_iso_past_cap <- percent(length(which(results_all_diagnosed_likelihood$maxDD>isolation_capacity_count))/length(results_all_diagnosed_likelihood$maxDD))
    
    results_all_diagnosed4 <- results_all_diagnosed3 %>%
      group_by(time2) %>%
      summarise_at(vars(maxD), list(nmin=min, Q1=~quantile(., probs = 0.25), Q95l=~quantile(., probs = 0.05),
                                    median=median, Q3=~quantile(., probs = 0.75),Q95u=~quantile(., probs = 0.95),
                                    max=max))
    
    

    avg_isolated_plot <- ggplot(data=results_all_diagnosed3, aes(x=time2, y=maxD, group=run)) + geom_line(alpha=0.3) +
      theme_classic() + theme(legend.position = "none") + 
      #geom_ribbon(aes(ymin = Q95l, ymax = Q95u), alpha = 0.1) + 
      xlab("Days") + ylab("Average number of isolated students") +
      ggtitle("Number of isolated \nstudents by day")
    
    
    # graph of presymptomatic over time
    
    presymp_plot <- ggplot(data=results_all_presymptomatic, aes(x=time, y=P_h, color=run)) + geom_line(alpha=0.3) +
        theme_classic() + theme(legend.position = "none") + ylab("Number infected presymptomatic") +
        xlab("Days")
    
    
    # graph of infected over time
    
    infected_plot <- ggplot(data=results_all_infected, aes(x=time, y=I_h, color=run)) + geom_line(alpha=0.3) +
        theme_classic() + theme(legend.position = "none") + ylab("Number infected symptomatic") +
        xlab("Days")
    
    
    # graph of new diagnoses
    
    ggplot(data=results_all_newlydiagnosed, aes(x=time, y=Dx0_h, color=run)) + geom_line(alpha=0.3)+
        theme_classic() + theme(legend.position = "none") + ylab("Number newly diagnosed") +
        xlab("Days")
    
    
    # graph of diagnosed (and isolated) over time
    
    diagnosed_plot <- ggplot(data=results_all_diagnosed, aes(x=time, y=Dx_h, color=run)) + geom_line(alpha=0.3)+
        theme_classic() + theme(legend.position = "none") + ylab("Number in isolation (diagnosed+)") +
        xlab("Days")
    
    
    #graph of recovered over time
    
    recovered_plot <- ggplot(data=results_all_recovered, aes(x=time, y=R_h, color=run)) + geom_line(alpha=0.3)+
        theme_classic() + theme(legend.position = "none") + ylab("Number recovered (cumulative)") +
        xlab("Days")
    
    
    # graph of all vaccinated over time
    
    vaccinated_plot <- ggplot(data=results_all_vaccinated, aes(x=time, y=Vs_h, color=run)) + geom_line(alpha=0.3)+
        theme_classic() + theme(legend.position = "none") + ylab("Number vaccinated") +
        xlab("Days")
    
    mean(results_all_vaccinated$Vs_h)
    min(results_all_vaccinated$Vs_h)
    max(results_all_vaccinated$Vs_h)
    
    
    output$DPlot <- renderPlot(avg_isolated_plot, height=300, width=300)
    output$QPlot <- renderPlot(avg_vaccinate_plot, height=300, width=300)
    output$IPlot <- renderPlot(avg_infectious_plot, height=300, width=300)
    output$D1Plot <- renderPlot(diagnosed_plot, height=300, width=300)
    output$Q1Plot <- renderPlot(vaccinated_plot, height=300, width=300)
    output$I1Plot <- renderPlot(infected_plot, height=300, width=300)
    output$R1Plot <- renderPlot(recovered_plot, height=300, width=300)
    # output$maxvaxplot <- renderPlot(maxvaxplot, height=300, width=300)
    # output$maxinfectionsplot <- renderPlot(maxinfectionsplot, height=300, width=300)
    # output$maxisoplot <- renderPlot(maxisoplot, height=300, width=300)
    #output$total_inf_hist <- renderPlot(total_inf_hist, height=200, width=200)
    #output$total_iso_hist
    
    #multiplot1(avg_isolated_plot, avg_vaccinate_plot, avg_infectious_plot,vaccinated_plot,diagnosed_plot,infected_plot,recovered_plot, cols=2)
    
    #multiplot2(avg_isolated_plot, avg_vaccinate_plot, avg_infectious_plot, cols=2)
    

    
    output$meaninfections_hr <- renderValueBox({
        valueBox(paste(round(additional_cases_mean_hr,0),"[",min_infections_hr,",",max_infections_hr,"]"), "Mean [min,max] number infections in 100 days, \ngiven initial case infects others, \nhigh-risk sub-group")
    })
    
    output$meaninfections_lr <- renderValueBox({
      valueBox(paste(round(additional_cases_mean_lr,0),"[",min_infections_lr,",",max_infections_lr,"]"), "Mean [min,max] number infections in 100 days, \nlow-risk sub-group")
    })
    
    output$meaninfections <- renderValueBox({
      valueBox(paste(round(cases_mean,0)), "Mean number infections in 100 days")
    })
    

    output$nonewinfections_hr <- renderValueBox({
        valueBox(additional_cases_hr, "Likelihood of additional infections following initial infections, \nin high-risk sub-group")
    })
    
    output$nonewinfections_lr <- renderValueBox({
      valueBox(additional_cases_lr, "Likelihood of infections in low-risk sub-group")
    })
    
    
    output$meanvaxcost <- renderValueBox({
      valueBox(mean(round(results_all_vaccinated_cost$cost,0)), "Mean cost of vaccination")
    })
    
    
    output$meandetectcost <- renderValueBox({
      valueBox(mean(round(results_all_diagnosed_cost3$cost,0)), "Mean cost of case isolation")
    })
    
    
    
})
    
     })
