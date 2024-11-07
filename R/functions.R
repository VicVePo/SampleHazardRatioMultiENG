#' Sample Size and Power Calculation for Survival Analysis (Cox Model)
#' 
#' @param alpha Significance level
#' @param power Desired power (optional if n is provided)
#' @param n Sample size (optional if power is provided)
#' @param HR Hazard Ratio to detect
#' @param CI_upper Upper confidence interval limit of HR
#' @param CI_lower Lower confidence interval limit of HR (optional if SE is provided)
#' @param SE Standard error of coefficient (optional if CI is provided)
#' @param n_previous Sample size of previous study (required for method "se")
#' @param sigma_x Standard deviation of predictor X (required for method "rho")
#' @param psi Probability that an observation is not censored (required for method "rho")
#' @param method Calculation method: "rho" or "se" (standard error)
#' @param rho_values Vector of rho values for sample size calculation (if method = "rho")
#' @return A data frame with results according to chosen method
#' @export
SampleHazardRatioMultiENG <- function(alpha, power = NULL, n = NULL, 
                                   HR = NULL,
                                   CI_upper = NULL, CI_lower = NULL,
                                   SE = NULL, n_previous = NULL,
                                   sigma_x = NULL, psi = NULL,
                                   method = "rho",
                                   rho_values = seq(0, 0.9, by = 0.1)) {
  
  # Verifications according to method
  if (method == "rho") {
    if (is.null(sigma_x) || is.null(psi)) {
      stop("For 'rho' method, sigma_x and psi are required")
    }
    if (is.null(HR)) {
      stop("HR must be specified")
    }
  } else if (method == "se") {
    if (is.null(HR)) {
      stop("For 'se' method, HR is required")
    }
    if (is.null(SE) && (is.null(CI_upper) || is.null(CI_lower))) {
      stop("For 'se' method, SE or both CI limits are required")
    }
    if (is.null(n_previous)) {
      stop("For 'se' method, n_previous is required")
    }
  } else {
    stop("Method must be 'rho' or 'se'")
  }
  
  # If CI provided, calculate SE
  if (!is.null(CI_upper) && !is.null(CI_lower)) {
    SE <- (log(CI_upper) - log(HR)) / qnorm(0.975)
    cat("Calculated Standard Error:", SE, "
")
  }
  
  if (method == "rho") {
    beta_a <- log(HR)
    z_alpha <- qnorm(1 - alpha/2)
    
    if (!is.null(power)) {
      z_gamma <- qnorm(power)
      
      calculate_n <- function(rho) {
        n <- ((z_alpha + z_gamma)^2) / ((beta_a * sigma_x)^2 * psi * (1 - rho^2))
        return(ceiling(n))
      }
      
      results <- data.frame(
        rho = rho_values,
        sample_size = sapply(rho_values, calculate_n)
      )
      
      results$expected_events <- ceiling(results$sample_size * psi)
    } else {
      calculate_power <- function(rho) {
        z_gamma <- sqrt(n * (beta_a * sigma_x)^2 * psi * (1 - rho^2)) - z_alpha
        power <- pnorm(z_gamma)
        return(power)
      }
      
      results <- data.frame(
        rho = rho_values,
        power = sapply(rho_values, calculate_power)
      )
      
      results$sample_size <- n
      results$expected_events <- ceiling(n * psi)
    }
  } else {
    # Standard error based method
    z_alpha <- qnorm(1 - alpha/2)
    if (!is.null(power)) {
      z_gamma <- qnorm(power)
      n <- ceiling((z_alpha + z_gamma)^2 * n_previous * SE^2 / (log(HR))^2)
      results <- data.frame(
        sample_size = n
      )
    } else {
      z_gamma <- sqrt(n * (log(HR))^2 / (n_previous * SE^2)) - z_alpha
      power <- pnorm(z_gamma)
      results <- data.frame(
        power = power,
        sample_size = n
      )
    }
  }
  
  return(results)
}

#' Logistics Calculation for Survival Studies
#'
#' @param total_n Total required sample size
#' @param expected_events Number of expected events
#' @param recruitment_rate Recruitment rate per month
#' @param followup_time Follow-up time in months
#' @param loss_to_followup_rate Expected loss to follow-up rate
#' @param rejection_rate Expected rejection rate
#' @return A list with study logistics calculations
#' @export
survival_study_logistics <- function(total_n, 
                                   expected_events,
                                   recruitment_rate,
                                   followup_time,
                                   loss_to_followup_rate,
                                   rejection_rate) {
  
  # Adjustment for losses to follow-up
  sample_with_losses <- total_n / (1 - loss_to_followup_rate)
  
  # Adjustment for rejections
  sample_to_contact <- sample_with_losses / (1 - rejection_rate)
  
  # Time calculations
  recruitment_time <- ceiling(sample_to_contact / recruitment_rate)
  total_study_time <- recruitment_time + followup_time
  
  # Results
  results <- list(
    sample_size = total_n,
    expected_events = expected_events,
    sample_with_losses = ceiling(sample_with_losses),
    sample_to_contact = ceiling(sample_to_contact),
    recruitment_months = recruitment_time,
    total_study_months = total_study_time
  )
  
  # Print summary
  cat("
Study logistics summary:
")
  cat("----------------------------------------
")
  cat("Required sample size:", total_n, "
")
  cat("Expected events:", expected_events, "
")
  cat("Sample size accounting for losses:", ceiling(sample_with_losses),
      "(", loss_to_followup_rate*100, "% losses)
")
  cat("Sample size to contact:", ceiling(sample_to_contact),
      "(", rejection_rate*100, "% rejection)
")
  cat("Recruitment months:", recruitment_time,
      "(", recruitment_rate, "persons per month)
")
  cat("Total study months:", total_study_time,
      "(including", followup_time, "months of follow-up)
")
  cat("----------------------------------------
")
  
  return(invisible(results))
}
