library(adaptivetau)
library(deSolve)
library(tidyverse)
library(ggnewscale)
library(future)
library(future.apply)
library(R.utils)

# resources ODE solving
RW_ode <- function(W0, r, D, media, k, tau, alpha = 1) {
  # Define the ODE
  ode_func <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dW <- W * monod(R, r, k)
      dR <- -dW * alpha
      return(list(c(dW = dW, dR = dR)))
    })
  }

  # Solve the ODE
  parameters <- c(r = r, D = D, media = media, k = k, tau = tau, alpha = alpha)
  out <- ode(y = c(W = W0, R = media - W0), times = c(0, tau),
    func = ode_func, parms = parameters)
  return(out[nrow(out), "W"])
}

# Find the equilibrium final bacterial population
find_W <- Vectorize(function(r, D, media, k, tau, flow = 1, alpha = 1) {
  # Define the function to be passed to uniroot
  func <- function(W0) log(W0 / RW_ode(W0, r, D, media, k, tau, alpha) / D)
  if (tau == 0 | D == 1) return(media - k/(r / flow - 1))
  if (tau < -log(D) / r * (1 + k / media)) return(0)
  # Use uniroot to find the root
  root <- uniroot(func, interval = c(1, media - 1))
  return(root$root / D)
}
)

# Probability a new mutant at the beginning of a growth phase will go extinct
phi <- function(D, s) {
    return(ifelse(D == 1, 1 / (1 + s), (1 - D) / (D ^ -s - D)))
}

# a growth rate function for nutrient-limited growth
monod <- function(R, r, k) {
  return(r * ifelse(k == 0, 1, 1 / (1 + k / R)))
}

# fitness advantages of new mutants
mutant_fitness <- function(num_mutants, r, w, names) {
  return(setNames(r * (1 + c(0, rexp(num_mutants, 1 / w))), names))
}

# a function that reduces all populations by a factor of D, in expectation
bottleneck <- function(state, config) {
  with(config, {
    pops <- setNames(rbinom(length(names), state[names], D), names)
    R <- state["R"] * D + media * (1 - D)
    return(c(pops, R))
  })
}

# a function outputting the transitions that can occur in the model
make_transitions <- function(names) {
    # Create the initial list
    rate_list <- list()
    # Add the growth rates for each mutant
    for (i in names) {
      rate_list[[paste0(i, "_growth")]] <- setNames(c(+1), i)
    }
    # Add the nutrient inflow rate
    rate_list[["R_inflow"]] <- setNames(c(+1), "R")
    # Add the death rates for each mutant
    for (i in names) {
      rate_list[[paste0(i, "_death")]] <- setNames(c(-1), i)
    }
    # Add the nutrient depletion rate
    rate_list[["R_depletion"]] <- setNames(c(-1), "R")
  return(rate_list)
}

# compute the rate at which each transition occurs
rates <- function(state, config, t) {
  with(as.list(c(state, config)), {
    # Calculate replication rates
    replication_rates <- state[names] * monod(R, r_vec, k)
    if (is.numeric(num_mutants)) {
        # find which mutant category we should be filling
        to_mutate <- min(which(state == 0), num_mutants + 2) - 1
        # chance of a replication in row i resulting in a strain j cell
        mutation <- diag(num_mutants + 1)
        if (to_mutate <= num_mutants) {
          mutation[1, to_mutate + 1] <- mu
          mutation[1, 1] <- 1 - mu
        } else {
          print("Out of mutant spots! :(")
        }
    }
    # Calculate growth rates including mutations
    growth_rates <- replication_rates %*% mutation
    # Calculate the outflow rates
    outflow_rates <- flow * state[names]
    # Calculate nutrient depletion rate
    R_depletion <- sum(replication_rates * alpha) + flow * R
    # Combine all rates and return
    rate_list <- c(growth_rates, flow * media, outflow_rates, R_depletion)
    return(setNames(rate_list, rep(c(names, "R"), 2)))
  })
}

# a function to implement one run of the model
single_run <- function(config, x) {
  with(config, {
    # Define the transitions of the model
    transitions <- make_transitions(names)
    # Initialise the state variables
    state <- init
    if (is.numeric(num_mutants)) {
      config$r_vec <- mutant_fitness(num_mutants, r, w, names)
    } else {
      config$r_vec <- r * t(1 + as.matrix(genotypes) %*% rexp(loci, 1 / w))
    }
    bottlenecks <- unique(round(c(seq(0, time, tau), time), 10))
    for (t in bottlenecks[-length(bottlenecks)]) {
      # Determine the time until the next bottleneck or dose
      end <- min(bottlenecks[bottlenecks > t] - t)
      # Create a new vector of growth rates for mutants, if locus-agnostic
      if (is.numeric(num_mutants)) {
        config$r_vec <- ifelse(state[0:num_mutants + 1] == 0,
          mutant_fitness(num_mutants, r, w, names), config$r_vec)
      }
      # set the seed for reproducibility
      if (is.numeric(seed)) set.seed(round(seed + (x * time + t) / tau))
      # Run the model between bottlenecks
      new <- ssa.adaptivetau(
        state, transitions, rates, config, tf = end,
        tl.params = list(maxtau = max_step),
        deterministic = grep("depletion", names(transitions))
      )
      # Make the time column reflect the overall time accurately
      new[, "time"] <- new[, "time"] + t
      # Avoid duplicate times by adding a small increment after the bottleneck
      new[1, "time"] <- new[1, "time"] * (1 + 1e-6)
      # Update the solution
      solution <- if (t == 0) new else rbind(solution, new)
      # Run the bottleneck and update the state
      state <- bottleneck(new[nrow(new), ], config)
    }
    # Interpolate the solution to the common time grid
    approx_vars <- lapply(colnames(solution), function(var) {
      approx(solution[, "time"], solution[, var], xout = time_grid)$y
    })
    solution_interpolated <- data.frame(
      setNames(approx_vars, colnames(solution)),
      rep = x
    )
    if (summarize) {
      # count the number of each mutant at the endpoint
      long <- pivot_longer(solution_interpolated, cols = -c(time, rep), names_to = "variable")
      final_counts <- long %>%
          group_by(rep, variable) %>%
          filter(time == endpoint, !(variable %in% c("W", "R"))) %>%
          summarise(final_value = value, .groups = "keep") %>%
          mutate(p_fix = 1 - current_phi ^ final_value)
      return(final_counts)
    } else {
      # add a new row at the end with the r_vec
      solution_interpolated[nrow(solution_interpolated) + 1, ] <-
        c(pi * 1e6, (config$r_vec / r) - 1, 0, x)
      return(solution_interpolated)
    }
  })
}

# a function to simulate the model
simulate <- function(
  rep = 1, # number of runs of the simulation
  seed = NULL, # seed for reproducibility
  time = 100, # time to simulate, in hours
  dt = 0.1, # time step, in hours
  max_step = Inf, # SSA max step parameter
  tau = 3, # frequency of bottlenecks, in hours
  D = 0.1, # dilution ratio at bottlenecks
  media = 1e9, # media nutrient concentration
  mu = 1e-9, # rate of beneficial mutations per replication
  N = 1e9, # (1/D) of the initial wild type population
  num_mutants = NULL, # number of mutants
  loci = NULL, # number of mutable loci in the genome
  r = 1, # wild type growth rate with infinite resources
  w = 0.1, # mean fitness effect size of beneficial mutation
  k = 1e8, # [R] at half-max growth rate
  alpha = 1, # nutrients used per replication
  flow = 0, # chemostat flow rate
  equilibrate = TRUE, # whether to start the population at equilibrium
  summarize = FALSE # whether to return a summary of the final population
  ) {
  # Define the parameters of the model
  if (is.numeric(num_mutants)) {
    names <- c("W", paste0("M", 1:num_mutants))
  } else {
      if (is.null(loci)) stop("Must specify num_mutants XOR loci")
      names <- paste0("G", intToBin(1:2^loci - 1))
      # Generate all possible binary strings of length loci
      genotypes <- expand.grid(replicate(loci, c(0, 1), simplify = FALSE))
      genotypes <- setNames(rev(genotypes), paste0("M", 1:loci))
      # Initialize the mutation matrix
      mutation <- diag(1 - mu, 2 ^ loci)
      # Loop over all genotypes
      for (i in 1:2^loci) {
        for (j in 1:2^loci) {
          # If the Hamming distance is 1, set the mutation rate to be nonzero
          if (sum(genotypes[i, ] != genotypes[j, ]) == 1) {
            mutation[i, j] <- mu / loci
          }
        }
      }
  }
  if (equilibrate & alpha != 0 & k != 0) {
    # Calculate the equilibrium population size
    N <- find_W(r, D, media, k, tau, flow)
  }
  init <- setNames(c(round(N * D), rep(0, length(names) - 1), media - N * D), c(names, "R"))
  # find the time just after the last bottleneck
  time_grid <- seq(0, time, by = dt) # a common time grid for all runs
  endpoint <- ceiling((time - time %% tau) / dt) * dt
  # find the likelihood of a new mutation at t=0 going extinct
  current_phi <- phi(D, w)
  config <- as.list(environment())
  # Run the simulation rep number of times, using parallelisation if possible
  plan(multisession) # compatible with both unix and Windows
  set.seed(seed) # set the seed for reproducibility
  long <- bind_rows(future_lapply(1:rep, function(x) {
    single_run(config, x)},
    future.seed = TRUE))
  if (!summarize) {
    # Convert the solutions to long format
    long <- pivot_longer(long, cols = -c(time, rep), names_to = "variable")
    config$s_all <- long %>% filter(near(time, pi * 1e6)) %>% select(-time)
    long <- long[!near(long$time, pi * 1e6), ]
  }
  return(list(long, config))
}

# Generate and analyse data
run_sims <- function(summary, rep = 1, time = 50, w = 0.1, r = 1, mu = 1e-9, dt = 1e-2,
    res = TRUE, flow = 1, num_mutants = 1e2, loci = NULL, summarize = TRUE) {
    for (i in seq_len(nrow(summary))) {
        D <- summary$D[i]
        tau <- summary$tau[i]
        k_ratio <- res *
            ifelse("k_ratio" %in% names(summary), summary$k_ratio[i], 1)
        data <- simulate(
            seed = i,
            rep = rep,
            time = time,
            dt = dt,
            tau = ifelse(tau == 0, time, tau),
            max_step = Inf,
            D = D,
            flow = ifelse(D == 1, flow, 0),
            media = 1e9,
            k = 1e9 * k_ratio,
            alpha = 1 * res,
            r = r * (1 + k_ratio),
            w = w,
            mu = mu,
            N = 1e9,
            num_mutants = num_mutants,
            loci = loci,
            summarize = summarize
        )
        if (summarize) {
            fixed <- data[[1]] %>%
                group_by(rep) %>%
                summarise(n = sum(final_value > 1e1), n_hat = sum(p_fix))
            # estimate the fixation rate and store this
            fixation_rate <- fixed$n_hat / data[[2]]$endpoint
            se <- sd(fixation_rate) / sqrt(length(fixation_rate))
            ci <- mean(fixation_rate) + se * qnorm(c(0.5, 0.025, 0.975))
        } else {
            func <- ifelse(is.null(loci), metric, metric_ci)
            ci <- func(data)
        }
        summary[i, c("rate", "ci_lower", "ci_upper")] <- ci
        print(i / nrow(summary))
    }
    return(summary)
}

# theoretical resource constrained adaptation rate
theory <- function(D, tau, r, w = 0.1, media = 1e9, k = 1e9, mu = 1e-9, flow = 1, alpha = 1) {
    k <- rep(k, length(D))[seq_along(D)]
    N <- ifelse(k == 0 | alpha == 0, media, find_W(r, D, media, k, tau, flow, alpha))
    x <- ifelse(D == 1 | tau == 0, flow, log(D) ^ 2 / (1 / D - 1) / tau)
    return(w / (1 + w) * mu * N * x)
}

# rate at a given time within [0, tau)
rate_at_t <- function(D, r, w, t) {
    u <- D * exp(r * t)
    theta <- (1 / D - 1) / (1 - D ^ w)
    return(r^2 / -log(D) * (u / (1 + theta * u ^ (1 + w))))
}

metric_ci <- function(data) {
    # find the number of each genotype at the endpoint
    final <- data[[1]] %>%
        filter(time == max(time), variable != "R")
    # Note which mutations each genotype has
    for (i in 1:data[[2]]$loci) {
        final[[paste0("M", i)]] <- substring(final$variable, i + 1, i + 1) == "1"
    }
    # Calculate the abundance and probability for each mutation
    counts <- final %>%
        pivot_longer(
            cols = starts_with("M"),
            names_to = "Mutant",
            values_to = "Mutant_value"
        ) %>%
        group_by(rep, Mutant) %>%
        summarise(p = sum(value * Mutant_value) / sum(value), .groups = "keep")
    # estimate the fixation rate and store this
    se <- sd(counts$p) / sqrt(length(counts$p))
    mean <- mean(counts$p)
    vec <- mean + se * qnorm(c(0.5, 0.025, 0.975))
    return(vec)
}
