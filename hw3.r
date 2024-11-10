# Load data
data <- read.table("hcmv.txt", header = TRUE)

# Parameters
num_simulations <- 10
num_palindromes <- 296
sequence_length <- 229354

# Generate simulated data
simulated_data <- list()
for (i in 1:num_simulations) {
    simulated_data[[i]] <- sample(1:sequence_length, num_palindromes, replace = FALSE)
}

# Function to calculate distances
calculate_distances <- function(data) {
    sorted_data <- sort(data)
    diff(sorted_data)
}

# Real and simulated distances
sorted_data <- sort(data$location)
real_distances <- diff(sorted_data)
simulated_distances <- lapply(simulated_data, calculate_distances)

# Plot real vs. simulated data
plot(data$location, rep(1, num_palindromes), pch = 16, col = "blue", main = "Real vs. Simulated Palindrome Locations", xlab = "Position along DNA sequence", ylab = "Simulation Index", ylim = c(0.5, num_simulations + 1))
for (i in 1:num_simulations) {
    points(simulated_data[[i]], rep(i + 1, num_palindromes), pch = 16, col = "red")
}

# Question 2a
hist(real_distances, breaks = 20, col = "blue", xlab = "Spacing", xlim = c(min(real_distances), max(unlist(simulated_distances))), ylab = "Frequency", freq = TRUE)
hist(simulated_distances[[1]], breaks = 20, col = rgb(1, 0, 0, 0.5), xlab = "Spacing", xlim = c(min(real_distances), max(unlist(simulated_distances))), ylab = "Frequency", freq = TRUE)

# Question 2b
palindrome_pairs_sum <- rowSums(cbind(sorted_data[-length(sorted_data)], sorted_data[-1]))
hist(palindrome_pairs_sum, breaks = 20, col = "blue", xlab = "Sum of Consecutive Pairs", ylab = "Frequency")

# Question 2c
num_triplets <- floor(length(sorted_data) / 3) * 3
sorted_data <- sorted_data[1:num_triplets]
triplets <- cbind(
    sorted_data[seq(1, num_triplets - 2, by = 3)],
    sorted_data[seq(2, num_triplets - 1, by = 3)],
    sorted_data[seq(3, num_triplets, by = 3)]
)
palindrome_triplets_sum <- rowSums(triplets)

# Simulate random palindrome locations
set.seed(123)
num_simulations <- 1000
random_simulations <- list()
for (i in 1:num_simulations) {
    simulated_locations <- sample(1:sequence_length, num_palindromes, replace = FALSE)
    sorted_simulated_locations <- sort(simulated_locations)
    num_sim_triplets <- floor(length(sorted_simulated_locations) / 3) * 3
    sorted_simulated_locations <- sorted_simulated_locations[1:num_sim_triplets]
    simulated_triplets <- cbind(
        sorted_simulated_locations[seq(1, num_sim_triplets - 2, by = 3)],
        sorted_simulated_locations[seq(2, num_sim_triplets - 1, by = 3)],
        sorted_simulated_locations[seq(3, num_sim_triplets, by = 3)]
    )
    simulated_triplets_sum <- rowSums(simulated_triplets)
    random_simulations[[i]] <- simulated_triplets_sum
}

# Plot observed and simulated triplets sums
hist(palindrome_triplets_sum, breaks = 30, col = rgb(0, 0, 1, 0.5), xlab = "Sum of Consecutive Triplets", ylab = "Frequency")
par(mfrow = c(2, 2))
for (i in 1:4) {
    hist(random_simulations[[i]], breaks = 30, col = rgb(1, 0, 0, 0.5), xlab = "Sum of Consecutive Triplets", ylab = "Frequency")
}
par(mfrow = c(1, 1))

# Question 3
interval_lengths <- c(1000, 10000)
count_palindromes <- function(positions, interval_length) {
    num_intervals <- ceiling(sequence_length / interval_length)
    interval_counts <- integer(num_intervals)
    for (pos in positions) {
        interval_index <- ceiling(pos / interval_length)
        interval_counts[interval_index] <- interval_counts[interval_index] + 1
    }
    return(interval_counts)
}

for (interval_length in interval_lengths) {
    par(mfrow = c(1, 2))
    num_intervals <- ceiling(sequence_length / interval_length)
    real_counts <- count_palindromes(data$location, interval_length)
    chisq_real <- chisq.test(real_counts, p = rep(1 / num_intervals, num_intervals))

    simulated_counts <- matrix(0, nrow = num_simulations, ncol = num_intervals)
    for (i in 1:num_simulations) {
        simulated_positions <- sample(1:sequence_length, num_palindromes, replace = FALSE)
        simulated_counts[i, ] <- count_palindromes(simulated_positions, interval_length)
    }

    hist(real_counts, breaks = 20, col = "blue", xlab = "Number of Palindromes", ylab = "Frequency")
    hist(simulated_counts[1, ], breaks = 20, col = rgb(1, 0, 0, 0.5), xlab = "Number of Palindromes", ylab = "Frequency")
}
par(mfrow = c(1, 1))

# Question 4
par(mfrow = c(1, 2))
interval_length <- 1000
real_counts <- count_palindromes(data$location, interval_length)
max_count_real <- max(real_counts)
max_interval_real <- which(real_counts == max_count_real)

simulated_max_counts <- numeric(num_simulations)
for (i in 1:num_simulations) {
    simulated_positions <- sample(1:sequence_length, num_palindromes, replace = FALSE)
    simulated_counts <- count_palindromes(simulated_positions, interval_length)
    simulated_max_counts[i] <- max(simulated_counts)
}

p_value <- mean(simulated_max_counts >= max_count_real)

plot(real_counts, type = "h", col = ifelse(real_counts == max_count_real, "red", "blue"), xlab = "Interval", ylab = "Palindrome Count")
hist(simulated_max_counts, breaks = 20, col = "gray", xlab = "Max Count", ylab = "Frequency")
abline(v = max_count_real, col = "red", lwd = 2, lty = 2)
par(mfrow = c(1, 1))

# Further Research Question
set.seed(42)
palindrome_lengths <- sample(10:18, num_palindromes, replace = TRUE)
interval_length <- 1000
real_counts <- count_palindromes(data$location, interval_length)
threshold <- median(real_counts)
high_density_intervals <- which(real_counts > threshold)
low_density_intervals <- which(real_counts <= threshold)

data$interval <- ceiling(data$location / interval_length)
high_density_lengths <- palindrome_lengths[data$interval %in% high_density_intervals]
low_density_lengths <- palindrome_lengths[data$interval %in% low_density_intervals]

t_test_result <- t.test(high_density_lengths, low_density_lengths, alternative = "two.sided")

par(mfrow = c(1, 2))
hist(high_density_lengths, breaks = 10, col = "blue", xlab = "Palindrome Length", ylab = "Frequency")
hist(low_density_lengths, breaks = 10, col = "red", xlab = "Palindrome Length", ylab = "Frequency")
par(mfrow = c(1, 1))
