# Function GET_FITTEDSINGLEDRUGS()
get_fittedsingledrugs <- function(drug_response, drug_name){	
	# Estimate the starting parameters ic50
    drug_response = as.data.frame(apply(drug_response, 2, as.numeric))
	estimate_param <- tryCatch({ 
        drm(effect ~ conc,
            data = drug_response,
            fct = LL.4( ## Use 4-parameter log-logistic model
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
                ##          HS    E_min E_inf  EC50
                ),
            na.action = na.omit,

            # uncomment these 2 lines to add boundaries to parameters
            #lowerl = c(0, 0, 0, 1e-6),
            #upperl = c(4, 100, 100, 1e+6),

            ## controlling constrained optimisation
            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
      }, warning = function(w){
            drm(effect ~ conc,
                data = drug_response,
                fct = LL.4( ## changed to log-logistic model
                    fixed = c(NA, NA, NA, NA),
                    names = c("SLOPE","MIN","MAX","IC50")
                    ##          HS    E_min E_inf  EC50
                    ),
                ## Reason for changing L.4 to LL.4:
                ##     Unless taking natural log(conc) when fitting,
                ##     otherwise the right form of the Hill model cannot be recovered.
                ##     Using L.4 is equivalent to fitting effect ~ exp(conc),
                ##     and the resulted EC50 is in fact exp(EC50)

                #logDose = 10, ## log based dose not used
                ## remove parameter boundaries due to their negative response values
                ## uncomment these 2 lines to add boundaries to parameters
                #lowerl = c(0, 0, 0, 1e-6),
                #upperl = c(4, 1, 1, 1e+6),

                na.action = na.omit)
      }, error = function(e){
            drm(effect ~ conc,
                data = drug_response,
                fct = LL.4( ## changed to log-logistic model
                    fixed = c(NA,NA,NA,NA),
                    names = c("SLOPE","MIN","MAX","IC50")
                    ##          HS    E_min E_inf  EC50
                ),
                na.action = na.omit,

                ## remove parameter boundaries due to their negative response values

                #logDose = 10 ## no log base 10 conc used
                control = drmc( ## set control argument
                    errorm = FALSE ## no error or warning on failed convergence
                )
            )
    })

    if (!is.null(estimate_param$convergence)) {
        message(drug_name, "failed to converge in fitting single-agent Hill curve.")
        return(NULL)
    } else {
        return(fitted(estimate_param))## like pred(), get predicted response from cruve
    }
}
c.get_fittedsingledrugs <- cmpfun(get_fittedsingledrugs)

# Function SDBASELINECOR() 
# for baseline correction using single drugs
# What is baseline correction???
SDbaselineCor <- function(plate.mat, conc.range, pair.list){
    conc.range <- apply(conc.range, 2, as.character)
    pair.list$index <- as.character(pair.list$index)
    pm <- as.matrix(plate.mat)
    pm <- apply(pm, 2, as.numeric)
    pm[nrow(plate.mat), ncol(plate.mat)] <- NA ## ???Why
    
    drugpair <- subset(pair.list, pair.list$index == 1) ##???
    #print(drugpair)
    d1c <- subset(
        conc.range,
        conc.range[,1] == drugpair$drug1
    )[, 2:ncol(conc.range)]
    #conc test range for drug 2 (2nd col in vv/pairslist)
    d2c <- subset(
        conc.range,
        conc.range[,1] == drugpair$drug2
    )[, 2:ncol(conc.range)] 

    ## create a viability matrix by doses
    rownames(pm) <- as.character(d2c)
    colnames(pm) <- as.character(d1c)
    
    # print(pm)
    d1 <- as.vector(unlist(rep(d1c,times = ncol(pm)))) 
    d2 <- as.vector(unlist(rep(d2c,times = ncol(pm)))) 
    
    ## repeat 6 concentrations by 6 times for each drug
    concd1 <- as.vector(sapply(d1c, rep, times = nrow(pm)))
    concd2 <- as.vector(sapply(d2c, rep, times = nrow(pm)))

    ## Effect of adding row drug to column drug
    coleffect <- as.vector(as.matrix(pm))
    ## Effect of adding column drug to row drug
    roweffect <- as.vector(t(as.matrix(pm))) 
    
    rowwise <- data.frame(as.numeric(d1), as.numeric(concd2), roweffect)
    # data combined row wise, rows are concatenated
    colwise <- data.frame(as.numeric(d2), as.numeric(concd1), coleffect)
    # data combined column wise, columns are concatenated
    
    ## Add drug 2 to drug 1
    colnames(rowwise)[c(1,2)] <- c("conc_d1", "conc_d2")
    ## Add drug 1 to drug 2
    colnames(colwise)[c(1,2)] <- c("conc_d2", "conc_d1")
    
    ## monotherapy response of drug 1
    single1 <- rowwise[1:nrow(plate.mat), c(1, 3)]
    # single drug denoted as the first row in the original matrix
    colnames(single1) <- c("conc","effect")
    # print(single1)
    
    ## monotherapy response of drug 2
    single2 <- colwise[1:ncol(plate.mat), c(1, 3)]
    # single drug denoted as the first column in the original matrix
    colnames(single2) <- c("conc","effect")
    
    # fitting single drugs independently
    ## dose-response curve fitting for drug 1
    sd1 <- get_fittedsingledrugs(
        drug_response = single1[-1, ], ## drop row 1 where drug 1 dose = 0
        drug_name = pair.list["drug1"]
    )
    #print(sd1)
    #print(single1)
    ## dose-response curve fitting for drug 2
    sd2 <- get_fittedsingledrugs(
        drug_response = single2[-1, ], ## drop row 1 where drug 2 dose = 0
        drug_name = pair.list["drug2"]
    )
    #print(sd2)
    #print(single2)
    
    if (is.null(sd1) & is.null(sd2)) {
        print(pair.list)
        return(NULL)
    }

    ## baseline = average of min predicted monotherapeutic responses
    if (!is.null(sd1) & is.null(sd2)) {
        message(paste("Drug 2", pair.list["drug2"], "failed to converge for obtaining baseline. Use drug1 only."))
        baseline <- min(as.numeric(sd1))
    } else if (is.null(sd1) & !is.null(sd2)) {
        message(paste("Drug 1", pair.list["drug1"], "failed to converge for obtaining baseline. Use drug2 only."))
        baseline <- min(as.numeric(sd2))
    } else {
        baseline <- (min(as.numeric(sd1)) + min(as.numeric(sd2))) / 2
    }

    
    ### Why did they do this
    # CORRECT MATRIX BY A WEIGHTED CORRECTION FACTOR
    pm_cor <- pm - ( (100 - pm) / 100 * baseline )
    
    rownames(pm_cor) <- as.character(d2c)
    colnames(pm_cor) <- as.character(d1c)
    
    output <- list(pm, pm_cor, drugpair)
    return(output) # return the adjusted matrix, and drug names
}
c.SDbaselineCor <- cmpfun(SDbaselineCor)

#####################################TWO WAY FITTING ###########################################################
# ----------------------------
# new functions defined here
# ----------------------------
twowayfitting <- function(cor_matrix, drug_pair, test_drc = FALSE) {
    #print(cor_matrix)
    #print(drug_pair)
    
    # Fitting single drugs using logistic functions
    # NA values treated
    # Drug 1 fitting (fitting the first row)
    #drug1=as.data.frame(mat.or.vec(7,0)) # first row in the matrix

    drug1 <- as.data.frame(
        mat.or.vec(
            nrow(cor_matrix) - 1,
            0
        )
    ) # first row in the matrix

    drug1$dose <- as.numeric(colnames(cor_matrix)[-1])
    #print(drug1)
    drug1$logconc <- log10(drug1$dose)
    drug1$inhibition <- as.numeric(cor_matrix[1, -1])

    if (test_drc) {
        ## For L.4
        drug1$lnconc <- log(drug1$dose)
        #drug1$lnlogconc <- log(drug1$logconc) logDose=10 => 10^(log(log10(dose))), non-sense
    }

    # NA values now treated

    drug1_model <- tryCatch({
        drm(
            #inhibition ~ logconc,
            inhibition ~ dose,
            data = drug1,
            fct = LL.4( ## changed to log-logistic model
                ## did not fix E_min and E_max in fitting single agent
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            #logDose = 10,

            # uncomment these 2 lines to add boundaries to parameters
            lowerl = ifelse(test_drc, c(0, 0, 0, 1e-6),
                            c(-Inf, -Inf, -Inf, -Inf)), ## remove boundaries in reproduction
            upperl = ifelse(test_drc, c(4, 100, 100, 1e+6),
                            c(Inf, Inf, Inf, Inf)),
            na.action = na.omit,

            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    }, error = function(e) {
        drm(
            #inhibition ~ logconc,
            inhibition ~ dose,
            data = drug1,
            fct = LL.4( ## changed to log-logistic model
                ## did not fix E_min and E_max in fitting single agent
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            #logDose = 10,
            # remove boundary
            na.action = na.omit,
            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    })
    
    if (test_drc) {

    ## Fitting log10(conc) vs natural scale
    ## to make sure my understanding of using drm with log10 base dose is right
    drug1_model_log10 <- tryCatch({
        drm(
            inhibition ~ logconc,
            data = drug1,
            fct = LL.4( ## changed to log-logistic model
                ## did not fix E_min and E_max in fitting single agent
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            logDose = 10,

            # uncomment these 2 lines to add boundaries to parameters
            lowerl = c(0, 0, 0, 1e-6),
            upperl = c(4, 100, 100, 1e+6),
            na.action = na.omit,

            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    }, error = function(e) {
        drm(
            inhibition ~ logconc,
            data = drug1,
            fct = LL.4( ## changed to log-logistic model
                ## did not fix E_min and E_max in fitting single agent
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            logDose = 10,
            # remove boundary
            na.action = na.omit,
            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    })



    ## compare to result fitted using L4
    ## code copied from the original script
    drug1_model_L4 = drm(inhibition ~ logconc, data = drug1, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10,na.action=na.omit)
    ## This is fitting
    ## L.4(logconc) = E_inf + (
    ##      (E_min - E_inf) / (1 + exp{HS(10^logconc - EC50)})
    ## )
    ## What Hill curve should be:
    ## Hill(logconc) = E_inf + (
    ##      (E_min - E_inf) / (1 + exp{HS( ln(10^logconc)  - ln(EC50))})
    ## )
    ## which can only be done with LL.4. Taking ln(logconc) will result in
    ## 10^ln(logconc) rather than ln(10^logconc), which is still wrong.
    ## The results are not convertiable between LL.4 and L.4 when fitting with log10(dose).

    ## Fitting L.4 with ln(dose)
    ## compare to result fitted with natural log(conc) using L4
    drug1_model_L4_ln = drm(inhibition ~ lnconc, data = drug1, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),na.action=na.omit)

    drug1_model_L4_ln_bounded = drm(
        inhibition ~ lnconc,
        data = drug1,
        fct = L.4(fixed = c(NA, NA, NA,NA),
                  names = c("SLOPE","MIN","MAX","IC50")),
        na.action=na.omit,
        lowerl = c(0, 0, 0, 1e-6),
        upperl = c(4, 100, 100, 1e+6)
    )

    ## Compare L.4 result against unbounded result fitted using LL.4
    drug1_model_unbounded <- drm(
            inhibition ~ logconc,
            data = drug1,
            fct = LL.4( ## changed to log-logistic model
                ## did not fix E_min and E_max in fitting single agent
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            logDose = 10,
            na.action = na.omit,

            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    
        drug1_fitted_inhibition_L4 <- PR(drug1_model_L4, drug1$logconc)

        drug1_fitted_inhibition_log10 <- PR(drug1_model_log10, drug1$logconc)

    ## Fit quality
    Rsqr_LL4_bounded <- 1 - (var(drug1$inhibition - fitted(drug1_model))/var(drug1$inhibition))
    Rsqr_LL4_unbounded <- 1 - (var(drug1$inhibition - fitted(drug1_model_unbounded))/var(drug1$inhibition))
    Rsqr_L4 <- 1 - (var(drug1$inhibition - fitted(drug1_model_L4))/var(drug1$inhibition))


    ## == BlockID 466 drug 1 result ===========================================

    ## -- Fitting with parameter boundary -------------------------------------

    ## -- drug1_model: fitted with LL.4, normal dose scale --------------------

    #A 'drc' model.
    #
    #Call:
    #drm(formula = inhibition ~ dose, data = drug1, fct = LL.4(fixed = c(NA,     NA, NA, NA), names = c("SLOPE
    #", "MIN", "MAX", "IC50")), na.action = na.omit,     control = drmc(errorm = FALSE), lowerl = c(0, 0, 0, 1
    #e-06),     upperl = c(4, 100, 100, 1e+06))
    #
    #Coefficients:
    #SLOPE:(Intercept)    MIN:(Intercept)    MAX:(Intercept)   IC50:(Intercept)
    #             0.00              45.08             100.00              65.38

    ## -- drug1_model_log10: fitted with LL.4, log10 dose scale --------------------

    #A 'drc' model.
    #
    #Call:
    #drm(formula = inhibition ~ logconc, data = drug1, fct = LL.4(fixed = c(NA,     NA, NA, NA), names = c("SL
    #OPE", "MIN", "MAX", "IC50")), na.action = na.omit,     logDose = 10, control = drmc(errorm = FALSE), lowe
    #rl = c(0,         0, 0, 1e-06), upperl = c(4, 100, 100, 1e+06))
    #
    #Coefficients:
    #SLOPE:(Intercept)    MIN:(Intercept)    MAX:(Intercept)   IC50:(Intercept)
    #             0.00              45.08             100.00              65.38

    #### This shows that we've used LL.4 to fit log10 based dose correctly. ###

    ## -- drug1_model_L4_ln_bounded: fitted with L.4, natural log dose scale ---
    
    #A 'drc' model.
    #
    #Call:
    #drm(formula = inhibition ~ lnconc, data = drug1, fct = L.4(fixed = c(NA,     NA, NA, NA), names = c("SLOP
    #E", "MIN", "MAX", "IC50")), na.action = na.omit,     lowerl = c(0, 0, 0, 1e-06), upperl = c(4, 100, 100,
    #1e+06))
    #
    #Coefficients:
    #SLOPE:(Intercept)    MIN:(Intercept)    MAX:(Intercept)   IC50:(Intercept)
    #             0.00              45.74              99.33               4.18
    
    ## Notice that other parameter estimations are close except for IC50, this is because:
    ## > exp(4.18)
    ## 65.36585

    ## This is the approach SynergyFinder uses when LL.4 goes wrong.
    ## However users are not aware that the estimated IC50 is in natural log,
    ## neither do they warn users about this.

    ## > Rsqr_LL4_unbounded
    ## [1] 0
    ## The poor fit is due to the negative response values and the bounded estimations.

    ## -- Unbounded fitting results --------------------------------------------

    ## -- drug1_model_unbounded: fitted with LL.4, normal dose scale -----------

    #A 'drc' model.
    #
    #Call:
    #drm(formula = inhibition ~ logconc, data = drug1, fct = LL.4(fixed = c(NA,     NA, NA, NA), names = c("SL
    #OPE", "MIN", "MAX", "IC50")), na.action = na.omit,     logDose = 10, control = drmc(errorm = FALSE))
    #
    #Coefficients:
    #SLOPE:(Intercept)    MIN:(Intercept)    MAX:(Intercept)   IC50:(Intercept)
    #           -1.808            -16.029             98.848             23.378

    ## -- drug1_model_L4_ln: fitted with L.4, natural log dose scale ----------

    #A 'drc' model.
    #
    #Call:
    #drm(formula = inhibition ~ lnconc, data = drug1, fct = L.4(fixed = c(NA,     NA, NA, NA), names = c("SLOP
    #E", "MIN", "MAX", "IC50")), na.action = na.omit)
    #
    #Coefficients:
    #SLOPE:(Intercept)    MIN:(Intercept)    MAX:(Intercept)   IC50:(Intercept)
    #           -1.808            -16.035             98.848              3.152

    ## Notice that all the other estimations are close except for IC50, again:
    ## > exp(3.152)
    ## [1] 23.38278
    ## which is still close to the IC50 estimation by LL.4

    ## This is to show that we have a correct understanding of L.4 vs. LL.4

    ## -- drug1_model_L4: fitted with L.4, from the code with the paper --------

    #A 'drc' model.
    #
    #Call:
    #drm(formula = inhibition ~ logconc, data = drug1, fct = L.4(fixed = c(NA,     NA, NA, NA), names = c("SLO
    #PE", "MIN", "MAX", "IC50")), na.action = na.omit,     logDose = 10)
    #
    #Coefficients:
    #SLOPE:(Intercept)    MIN:(Intercept)    MAX:(Intercept)   IC50:(Intercept)
    #         -0.03976         -853.34093           97.67009          -45.80743

    ## This is the Boltzmann model, not the Hill model.
    ## Its IC50 cannot be transformed back to that of the Hill model
    ## because it was fitted using L.4 with x = 10^log10(dose) rather than x = ln(10^log10(dose))
    ## Its very negative E_min and IC50 are also uninterpretable.

    ## We have shown that the approach used in the original code with the paper is incorrect.
    } ## test_drc END

    ## No $convergence on successful convergence
    if (is.null(drug1_model$convergence)) {
       drug1$fitted_inhibition <- PR(drug1_model, drug1$dose)
    } else {
        message(paste("Drug 1", drug_pair["drug1"], "failed to converge in fitting single-agent dose response curve in 2-way fitting."))
        drug1$fitted_inhibition <- NA_real_
    }
    

    # Drug 2 fitting (fitting the first column)
    #drug2=as.data.frame(mat.or.vec(7,0)) # first column in the matrix
    drug2 <- as.data.frame(mat.or.vec(ncol(cor_matrix) - 1, 0))
    # first column in the matrix
    drug2$dose <- as.numeric(rownames(cor_matrix)[-1])
    #print(drug2)
    drug2$logconc <- log10(drug2$dose)
    drug2$inhibition=as.numeric(cor_matrix[-1,1])

    drug2_model <- tryCatch({
        drm(
            #inhibition ~ logconc,
            inhibition ~ dose,
            data = drug2,
            fct = LL.4( ## changed to log-logistic model
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            #logDose=10,

            # uncomment these 2 lines to add boundaries to parameters
            lowerl = ifelse(test_drc, c(0, 0, 0, 1e-6),
                            c(-Inf, -Inf, -Inf, -Inf)), ## remove boundaries in reproduction
            upperl = ifelse(test_drc, c(4, 100, 100, 1e+6),
                            c(-Inf, -Inf, -Inf, -Inf)),

            na.action=na.omit
        )
    }, error = function(e) {
        drm(
            #inhibition ~ logconc,
            inhibition ~ dose,
            data = drug2,
            fct = LL.4( ## changed to log-logistic model
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            #logDose=10,
            ## remove parameter boundary due to negative response values
            na.action=na.omit,
            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    })

    if (test_drc) {
    drug2_model_log10 <- tryCatch({
        drm(
            inhibition ~ logconc,
            data = drug2,
            fct = LL.4( ## changed to log-logistic model
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            logDose=10,

            # uncomment these 2 lines to add boundaries to parameters
            lowerl = c(0, 0, 0, 1e-6),
            upperl = c(4, 100, 100, 1e+6),

            na.action=na.omit
        )
    }, error = function(e) {
        drm(
            inhibition ~ logconc,
            data = drug2,
            fct = LL.4( ## changed to log-logistic model
                fixed = c(NA, NA, NA, NA),
                names = c("SLOPE", "MIN", "MAX", "IC50")
            ),
            logDose=10,
            ## remove parameter boundary due to negative response values
            na.action=na.omit,
            control = drmc( ## set control argument
                errorm = FALSE ## no error or warning on failed convergence
            )
        )
    })

    ## Compare against result fitted with L4:

    drug2_model_L4=drm(inhibition ~ log(dose), data = drug2, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),na.action=na.omit)

        drug2_fitted_inhibition_L4 <- PR(drug2_model_L4, drug2$logconc)
        drug2_fitted_inhibition_log10 <- PR(drug2_model_log10, drug2$logconc)
    }


    if (is.null(drug2_model$convergence)) {
       drug2$fitted_inhibition <- PR(drug2_model, drug2$dose)
    } else {
        message(paste("Drug 2", drug_pair["drug2"], "failed to converge in fitting single-agent dose response curve in 2-way fitting."))
        drug2$fitted_inhibition <- NA_real_
    }
    
    # Update the first row and first column
    cor_matrix_2 = mat.or.vec(nrow(cor_matrix), ncol(cor_matrix))
    colnames(cor_matrix_2) = colnames(cor_matrix)
    rownames(cor_matrix_2) = rownames(cor_matrix)
    cor_matrix_2[1, c(2:ncol(cor_matrix))] = drug1$fitted_inhibition
    cor_matrix_2[c(2:nrow(cor_matrix)), 1] = drug2$fitted_inhibition
    ## matrix containing only predicted response to single agents
    #print(cor_matrix_2)
    
    # Update the column2-column8
    ## Fitting response curve of adding column drug to row drug
    cor_matrix_3 = cor_matrix_2
    for (i in 2:ncol(cor_matrix)){
            tmp = as.data.frame(mat.or.vec(nrow(cor_matrix)-1,0))
            tmp$dose = drug2$dose
            tmp$logconc = drug2$logconc
            tmp$inhibition = cor_matrix[c(2:nrow(cor_matrix)), i]
            tmp_min = cor_matrix_2[1, i] ## Add column drug to row drug
            tmp_model <- tryCatch({
                drm(
                    #inhibition ~ logconc,
                    inhibition ~ dose,
                    data = tmp,
                    fct = LL.4( ## changed to log-logistic model
                        fixed = c(NA, tmp_min, NA,NA), ## Fix E_min = E_col(i)
                        names = c("SLOPE", "MIN", "MAX", "IC50")
                    ),
                    #logDose = 10,

                    # uncomment these 2 lines to add boundaries to parameters
                    #lowerl = c(0, 0, 1e-6),
                    #upperl = c(4, 100, 1e+6),

                    na.action = na.omit
                )
            }, error = function(e) {
                drm(
                    #inhibition ~ logconc,
                    inhibition ~ dose,
                    data = tmp,
                    fct = LL.4( ## changed to log-logistic model
                        fixed = c(NA, tmp_min, NA,NA), ## Fix E_min = E_col(i)
                        names = c("SLOPE", "MIN", "MAX", "IC50")
                    ),
                    #logDose = 10,
                    ## remove boundary due to their negative response
                    na.action = na.omit,
                    control = drmc( ## set control argument
                        errorm = FALSE ## no error or warning on failed convergence
                    )
                )
            })
            if (is.null(tmp_model$convergence)) {
                tmp$fitted_inhibition = PR(tmp_model, tmp$dose)
            } else {
                message(paste("Adding drug", drug_pair["drug1"], "failed to converge in 2-way fitting at dose", tmp_min))
                tmp$fitted_inhibition = NA_real_
            }
            if (!any(is.na(tmp$fitted_inhibition))){
                if(tmp$fitted_inhibition[ncol(cor_matrix)-1] < 0)
                    tmp$fitted_inhibition[ncol(cor_matrix)-1] <- tmp_min
                    ## Replace negative response with E_min(added drug)
            }
            cor_matrix_3[c(2:nrow(cor_matrix)),i] = tmp$fitted_inhibition
            #print(tmp)
    }
    #print(cor_matrix_3)
    
    # Update the row2-row8
    ## Fitting response curve of adding row drug to column drug
    cor_matrix_4 = cor_matrix_2
    for (i in 2:nrow(cor_matrix)){
    		tmp = as.data.frame(mat.or.vec(nrow(cor_matrix)-1,0))
    		tmp$dose = drug1$dose
    		tmp$logconc = drug1$logconc
    		tmp$inhibition = cor_matrix[i,c(2:ncol(cor_matrix))]
    		tmp_min = cor_matrix_2[i,1] ## Add row drug to column drug
    		tmp_model <- tryCatch({
                drm(
                    #inhibition ~ logconc,
                    inhibition ~ dose,
                    data = tmp,
                    fct = LL.4( ## changed to log-logistic model
                        fixed = c(NA, tmp_min, NA, NA), ## Fix E_min = E_row(i)
                        names = c("SLOPE", "MIN", "MAX", "IC50")
                    ),
                    #logDose = 10,

                    # uncomment these 2 lines to add boundaries to parameters
                    #lowerl = c(0, 0, 1e-6),
                    #upperl = c(4, 100, 1e+6),
                    
                    na.action = na.omit,
                    control = drmc( ## set control argument
                        errorm = FALSE ## no error or warning on failed convergence
                    )
                )
            }, warning = function(w) {
                drm(
                    #inhibition ~ logconc,
                    inhibition ~ dose,
                    data = tmp,
                    fct = LL.4( ## changed to log-logistic model
                        fixed = c(NA, tmp_min, NA, NA), ## Fix E_min = E_row(i)
                        names = c("SLOPE", "MIN", "MAX", "IC50")
                    ),
                    #logDose = 10,
                    ## remove boundary due to their negative response
                    na.action = na.omit,
                    control = drmc( ## set control argument
                        errorm = FALSE ## no error or warning on failed convergence
                    )
                )
            }, error = function(e) {
                drm(
                    #inhibition ~ logconc,
                    inhibition ~ dose,
                    data = tmp,
                    fct = LL.4( ## changed to log-logistic model
                        fixed = c(NA, tmp_min, NA, NA), ## Fix E_min = E_row(i)
                        names = c("SLOPE", "MIN", "MAX", "IC50")
                    ),
                    #logDose = 10,
                    ## remove boundary due to their negative response
                    na.action = na.omit,
                    control = drmc( ## set control argument
                        errorm = FALSE ## no error or warning on failed convergence
                    )
                )
            })
            if (is.null(tmp_model$convergence)) {
                tmp$fitted_inhibition = PR(tmp_model, tmp$dose)
            } else {
                message(paste("Adding drug", drug_pair["drug2"], "failed to converge in 2-way fitting at dose", tmp_min))
                tmp$fitted_inhibition = NA_real_
            }
            if (!any(is.na(tmp$fitted_inhibition))){
                if(tmp$fitted_inhibition[ncol(cor_matrix)-1] < 0)
                    tmp$fitted_inhibition[ncol(cor_matrix)-1] = tmp_min
                    ## Replace negative response with E_min(added drug)
            }
            cor_matrix_4[i,c(2:ncol(cor_matrix))] = tmp$fitted_inhibition
            #print(tmp)
     }
    #print(cor_matrix_4)
     
    ## See eqn. (19)
    # take average of cor_matrix_3 and cor_matrix_4 as cor_matrix_final
    cor_matrix_final = (cor_matrix_3+cor_matrix_4) / 2
    #print(cor_matrix_final)
    
    # make cor_matrix_bliss based on cor_matrix_2
    cor_matrix_bliss = cor_matrix_2
    for (i in 2:nrow(cor_matrix_2)){
    	for (j in 2:ncol((cor_matrix_2))){
    		cor_matrix_bliss[i,j] <- (
                cor_matrix_2[i,1] +
                cor_matrix_2[1,j] -
                cor_matrix_2[i,1] * cor_matrix_2[1,j] / 100
            )
    	}
    }
    # negative and positive controls are removed
    cor_matrix_final[1,1] <- 0
    cor_matrix_bliss[1,1] <- 0 
    #cor_matrix_bliss[8,8]=0
    cor_matrix_final[nrow(cor_matrix_final),
                     ncol(cor_matrix_final)] = ifelse(
                         cor_matrix_final[nrow(cor_matrix_final),
                                          ncol(cor_matrix_final)] > 100,
                         100,
                         cor_matrix_final[nrow(cor_matrix_final),
                                          ncol(cor_matrix_final)]
                     ) # cannot be over 100 for the estimation
                       ## So they truncated E_inf at 100
    
    ## See eqn. (19) for why taking diff
    diff_cor_matrix <- (
        cor_matrix_final ## Two way fitted result
        -
        cor_matrix_bliss ## Bliss product of predicted monotherapeutic responses
    )
    
    diff_cor_matrix_new <- diff_cor_matrix[-1,] ### remove first row
    diff_cor_matrix_new <- diff_cor_matrix[,-1] ### remove first col

    ## Average of delta scores over all dose combinations
    sum_diff_cor_matrix <- (
        sum(diff_cor_matrix_new, na.rm = TRUE) /
        ( ncol(diff_cor_matrix_new) * nrow(diff_cor_matrix_new) )
    )

    output = sum_diff_cor_matrix
    return(output)
}
#END of the twowayfitting function

c.twowayfitting <- cmpfun(twowayfitting)


