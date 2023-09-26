source("simulation.R")
library(scales)
library(scico)
library(gridExtra)

r <- 1 # resource unconstrained growth rate
r_adj <- 1.023 # adjustment for non-infinitesmal step size
r_res <- 1.5 # resource constrained growth rate
media <- 1e9 # resource concentration in dilution media
w <- 0.1 # average selective benefit
fig_dir <- "figs"
data_dir <- "data"

custom_theme <- theme(
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "bottom",
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    strip.text = element_text(size = 25),
    legend.key.size = unit(1.5, "cm"),
  )

base_plot <- function(summary) {
    plot <- ggplot(summary) +
        scale_x_continuous(trans = log10_trans(),
            breaks = 10^seq(-7, 0),
            labels = trans_format("log10", math_format(10^.x))) +
        scale_y_continuous(trans = log10_trans(),
            breaks = 10^seq(-7, 0),
            labels = trans_format("log10", math_format(10^.x))) +
        labs(
            x = expression(italic("D")),
            y = expression(paste("fixation rate (loci hour"^"-1",")"))
        ) +
        theme_light() +
        custom_theme
    return(plot)
}

dynamics <- function(data, part) {
    # Filter data to include only non-zero mutants and remove "R"
    non_zero_data <- data %>%
    group_by(variable) %>%
    filter(sum(value) > 0, variable != "R")

    # Define colors using scico
    unique_vars <- unique(non_zero_data$variable)
    colors <- setNames(scico(length(unique_vars), palette = "roma"), unique_vars)
    colors["W"] <- "black"

    # Create plot
    p <- ggplot(non_zero_data, aes(x = time, y = value, color = variable)) +
        geom_line(linewidth = 1.5) +
        scale_color_manual(values = colors) +
        scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
            breaks = 10^seq(0, 9),
            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(
            x = "time (hours)",
            y = "population size",
            color = "Strain",
            fill = "Strain"
        ) +
        annotate("text", x = 0, y = 3e9, size = 15,
            label = part, fontface = "bold", hjust = 0.1, vjust = 0.3) +
        theme_light() +
        custom_theme
    return(p)
}

constrained <- function(summary) {
    p <- ggplot(summary, aes(x = D, y = tau)) +
        geom_tile(aes(fill = log10(rate))) +
        scale_x_continuous(trans = scales::log10_trans(),
            breaks = 10^seq(-4, 0),
            labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        scale_y_continuous(trans = scales::log2_trans(),
            breaks = c(24, 6, 2, 0.5)) +
        scale_fill_gradient(low = "white", high = "blue",
            breaks = pretty_breaks(n = 2), labels = scales::math_format(10^.x)) +
        labs(x = expression(italic("D")),
            y = expression(paste(italic(tau), " (hours)")),
            fill = "fixation rate") +
        theme_minimal() +
        custom_theme
    return(p)
}

### fig dynamics
data <- list()
for (i in 1:5) {
    res <- i != 5
    k_ratio <- ifelse(res, 10 ^ (i - 3), 0)
    data[[i]] <- simulate(
    seed = 1,
    time = 30,
    tau = log(10),
    D = 0.1,
    media = 1e9,
    mu = 3e-9,
    N = 1e9,
    num_mutants = 9,
    r = r * r_adj * r_res ^ res * (k_ratio + 1),
    w = w,
    k = k_ratio * 1e9,
    alpha = res
)[[1]]
}

# save the plots
pdf(paste0(fig_dir, "/dynamics.pdf"), width = 10, height = 10)
dynamics(data[[5]], "")
dev.off()

pdf(paste0(fig_dir, "/dynamics-constrained.pdf"), width = 20, height = 20)
grid.arrange(
    dynamics(data[[1]], "A"),
    dynamics(data[[2]], "B"),
    dynamics(data[[3]], "C"),
    dynamics(data[[4]], "D"),
    ncol = 2
)
dev.off()

### fig:optimal
summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$tau <- - log(summary$D)
summary <- run_sims(summary, rep = 1e3, r = r * r_adj, w = w, res = FALSE)
summary$theory <- theory(summary$D, summary$tau, r, w, k = 0)

# Save the summary to a file
save(summary, file = paste0(data_dir, "/optimal.rdata"))
load(paste0(data_dir, "/optimal.rdata"))

# save the plot
pdf(paste0(fig_dir, "/optimal.pdf"), width = 10, height = 10)
base_plot(summary) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), linewidth = 0.6) +
    geom_point(aes(x = D, y = rate), size = 3) +
    geom_line(aes(x = D, y = theory), linewidth = 1)
dev.off()

### fig constrained

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), tau = 24 * 2 ^ - seq(0, 6.5, by = 0.5))
summary <- run_sims(summary, rep = 1e3, r = r_res, res = TRUE)

# Save the summary to a file
save(summary, file = paste0(data_dir, "/constrained.rdata"))
load(paste0(data_dir, "/constrained.rdata"))
# save the plot
pdf(paste0(fig_dir, "/constrained.pdf"), width = 10, height = 10)
constrained(summary)
dev.off()

theory_data <- summary
theory_data$rate <- theory(theory_data$D, theory_data$tau, r_res * 2)
comparison <- theory_data
comparison$rate <- theory_data$rate / summary$rate

pdf(paste0(fig_dir, "/constrained_comp.pdf"), width = 10, height = 20)
grid.arrange(
    constrained(summary) +
        labs(fill = "numerical") +
        annotate("text", x = 1e-4, y = 24, size = 15,
            label = "A", fontface = "bold", hjust = 0, vjust = 1),
    constrained(theory_data) +
        labs(fill = "analytical") +
        annotate("text", x = 1e-4, y = 24, size = 15,
            label = "B", fontface = "bold", hjust = 0, vjust = 1),
    constrained(comparison) +
        labs(fill = "analytical/numerical") +
        annotate("text", x = 1e-4, y = 24, size = 15,
            label = "C", fontface = "bold", hjust = 0, vjust = 1),
    ncol = 1
)
dev.off()

### fig:tau24

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1), tau = 24)
summary <- run_sims(summary, rep = 2e3, r = r_res, res = TRUE)
summary$theory <- theory(summary$D, summary$tau, r_res * 2)

# Save the summary to a file
save(summary, file = paste0(data_dir, "/tau24.rdata"))
load(paste0(data_dir, "/tau24.rdata"))
# save the plot
pdf(paste0(fig_dir, "/tau24.pdf"), width = 10, height = 10)
base_plot(summary) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), linewidth = 0.8) +
    geom_point(aes(x = D, y = rate), size = 3) +
    geom_line(aes(x = D, y = theory), linewidth = 1)
dev.off()

### fig:k_variation_optimal

summary <- expand.grid(D = 10 ^ - seq(0, 4, by = 0.1),
    k_ratio = 10 ^ seq(-2, 1, by = 1))
summary$tau <- - log(summary$D)
summary$k <- as.factor(summary$k_ratio * 1e9)
summary <- run_sims(summary, rep = 1e3, time = 20, r = r_res, res = TRUE)
summary$theory <- theory(summary$D, summary$tau, r_res * (1 + summary$k_ratio), k = 1e9 * summary$k_ratio)

# Save the summary to a file
save(summary, file = paste0(data_dir, "/k_variation_optimal.rdata"))
load(paste0(data_dir, "/k_variation_optimal.rdata"))
# save the plot
pdf(paste0(fig_dir, "/k_variation_optimal.pdf"), width = 10, height = 10)
base_plot(summary) +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper, color = k),
        linewidth = 0.8) +
    geom_line(aes(x = D, y = theory, color = k), linewidth = 1) +
    geom_point(aes(x = D, y = rate, color = k), size = 3) +
    scale_color_scico_d(palette = "roma")
dev.off()


### fig:t_distribution
D <- 10^-0.1
data <- simulate(seed = 1, time = 50, rep = 1e2, dt = 1e-2, max_step = 1e-2,
    D = D,  w = w, tau = -log(D), media = 1e9, N = 1e9,
    mu = 1e-8, k = 0, alpha = 0, r = r * r_adj, num_mutants = 2e2)
final_counts <- data[[1]] %>%
        group_by(rep, variable) %>%
        filter(time == data[[2]]$endpoint, !(variable %in% c("W", "R"))) %>%
        summarise(final_value = value, .groups = "keep") %>%
        mutate(p_fix = 1 - data[[2]]$current_phi ^ final_value)

mutation_times <- data[[1]] %>%
    filter(value == 0) %>%
    group_by(rep, variable) %>%
    summarise(last_zero = max(time)  %% data[[2]]$tau, .groups = "keep")
t_data <- final_counts %>%
    inner_join(mutation_times, by = c("rep", "variable"))

# save the data
save(t_data, file = paste0(data_dir, "/t_distribution.rdata"))

t_theory <- data.frame(
  t = seq(0, data[[2]]$tau, by = 0.001)
)
t_theory$rate <- rate_at_t(D, r = r, w = w, t = t_theory$t)
t_theory$rate <- t_theory$rate / sum(t_theory$rate * 0.001)

# save the plot
pdf(paste0(fig_dir, "/t_distribution.pdf"), width = 10, height = 10)
ggplot(t_data, aes(x = last_zero)) +
  geom_density(aes(y = after_stat(density), weight = p_fix),
    adjust = 1 / 2, fill = "grey", alpha = 1, bw = 0.01) +
  geom_line(data = t_theory, aes(x = t, y = rate, color = "theory")) +
  scale_color_manual(name = NULL, values = c("theory" = "blue")) +
  labs(
    x = expression(italic(t)),
    y = "probability density"
  ) +
  theme_minimal() +
  custom_theme
dev.off()

### fig:s_distribution
s_data <- final_counts %>%
    inner_join(data[[2]]$s_all, by = c("rep", "variable"))

# save the data
save(s_data, file = paste0(data_dir, "/s_distribution.rdata"))

s_theory <- data.frame(
  s = seq(0, 1, by = 0.001)
)
s_theory$fix <- s_theory$s * exp(-s_theory$s / w) / w^2
s_theory$arise <- exp(-s_theory$s / w) / w

# save the plot
pdf(paste0(fig_dir, "/s_distribution.pdf"), width = 10, height = 10)
ggplot(s_data, aes(x = value)) +
  # raise p_fix to a high power to ensure mutations actually very likely to fix
  # are the ones that dominate the distribution
  geom_density(aes(y = after_stat(density), weight = p_fix ^ 20),
    adjust = 1 / 2, fill = "grey", alpha = 1, bw = 0.01) +
  geom_line(data = s_theory, aes(x = s, y = fix, color = "fixing")) +
  geom_line(data = s_theory, aes(x = s, y = arise, color = "arising")) +
  scale_color_manual(name = "mutations",
    values = c("fixing" = "blue", "arising" = "red")) +
  labs(
    x = expression(italic(t)),
    y = "probability density"
  ) +
  theme_minimal() +
  custom_theme
dev.off()

# weighted mean of s for mutations that go onto fix
weighted.mean(s_data$value, s_data$p_fix ^ 20)

### fig:binomial
# define parameters
D <- 0.9
s <- 0.1
rep <- 1e6
label <- c("det_poi", "det_bin", "sto_poi", "sto_bin")
shapes <- 15:18
set.seed(0)

# define a tibble with the source and pmf values
pmfs <- tibble(
  source = rep(c("det_bin", "det_poi", "sto_bin", "sto_poi"), each = rep),
  value = c(
    {e <- D^-(1 + s)
    f <- floor(e)
    rbinom(rep, f, D) + rbinom(rep, 1, (e - f) * D)},
    rpois(rep, D^-s),
    {m <- 1 + rgeom(rep, D ^ (1 + s))
    rbinom(rep, m, D)},
    {m <- 1 + rgeom(rep, D ^ (1 + s))
    rpois(rep, m)}
  )
) %>%
  group_by(source, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(pmf = n / sum(n))
pmfs$source <- factor(pmfs$source, levels = c("det_poi", "sto_poi", "det_bin", "sto_bin"))

# save the plot
pdf(paste0(fig_dir, "/binomial.pdf"), width = 10, height = 10)
ggplot(pmfs, aes(x = value, y = pmf, colour = source, shape = source)) +
  geom_jitter(size = 5, width = 0.2, height = 0) +
  scale_color_manual(values = setNames(scico(4, palette = "roma"), label)) +
  scale_shape_manual(values = setNames(shapes, label)) +
  theme_light() +
  scale_x_continuous(breaks = 0:7, limits = c(-0.2, 7.2)) +
  scale_y_continuous(trans = scales::log10_trans(),
        breaks = 10^seq(-7, 0),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(
    x = expression(paste("mutants remaining after bottlenecking, ",
        italic(M(tau^"+")))),
    y = "probability mass",
    colour = NULL,
    shape = NULL
  ) +
  custom_theme
dev.off()


### fig:methodology

summary <- expand.grid(D = 10 ^ - seq(0.01, 4, 0.01), div_tau = c(FALSE, TRUE))
summary$tau <- -log(summary$D) / r
summary <- with(summary, {
    summary$det_bin <- D * (1 - phi(D, w)) * (tau ^ !div_tau)
    summary$det_poi <- 2 * w * D * (log(D)^2) / (tau ^ div_tau)
    summary$sto_bin <- theory(D, tau, r, w, k = 0) * (tau ^ !div_tau)
    summary$sto_poi <- D * (1 - D ^ w) * (tau ^ !div_tau)
    return(summary)
})

# pivot longer
summary <- summary %>%
    pivot_longer(cols = !c("D", "div_tau", "tau"),
        names_to = "model", values_to = "rate")
summary$model <- factor(summary$model, levels = c("det_poi", "sto_poi", "det_bin", "sto_bin"))
summary$div_tau <- factor(summary$div_tau, levels = c(TRUE, FALSE))
# save the plot
pdf(paste0(fig_dir, "/methodology.pdf"), width = 10, height = 10)
base_plot(summary) +
    geom_line(aes(x = D, y = rate, color = model, linetype = div_tau), linewidth = 1) +
    scale_color_manual(values = setNames(scico(4, palette = "roma"), label)) +
    scale_linetype_manual(values = c("solid", "dotted"), labels = c("time", "transfers")) +
    labs(color = "model", linetype = "optimize") +
    guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    theme(legend.text = element_text(size = 23),
        legend.title = element_text(size = 23))
dev.off()

### fig:ci

summary <- expand.grid(D = 10 ^ - seq(0.1, 4, by = 0.1))
summary$tau <- - log(summary$D)
summary <- run_sims(summary, rep = 2e2, time = 1e3, r = r_res,
    loci = 4, num_mutants = NULL, res = TRUE, dt = 10, summarize = FALSE)

# Save the summary to a file
save(summary, file = paste0(data_dir, "/ci.rdata"))

# save the plot
pdf(paste0(fig_dir, "/ci.pdf"), width = 10, height = 10)
base_plot(summary) +
    labs(y = "Average mutation frequency") +
    geom_errorbar(aes(x = D, ymin = ci_lower, ymax = ci_upper), linewidth = 0.8) +
    geom_point(aes(x = D, y = rate), size = 3)
dev.off()
