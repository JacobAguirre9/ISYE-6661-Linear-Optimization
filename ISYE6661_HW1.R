# HW 1 Script for ISYE 6661
# @JacobAguirre

#Plotting Q1

library(plotly)

# Define the constraints
con1 <- function(x1, x2) x1 - 2*x2
con2 <- function(x1, x2) 2*x1 - 3*x2 - 2
con3 <- function(x1, x2) x1 - x2 - 3
con4 <- function(x1, x2) -x1 + 2*x2 - 2
con5 <- function(x1, x2) -2*x1 + x2

# Create a grid of points
x1 <- seq(-4, 4, length.out = 100)
x2 <- seq(-4, 4, length.out = 100)
grid <- expand.grid(x1, x2)

# Check which points satisfy the constraints
# Check which points satisfy the constraints
df <- as.data.frame(grid)
df$satisfies_con1 <- con1(df$Var1, df$Var2) <= 0
df$satisfies_con2 <- con2(df$Var1, df$Var2) <= 0
df$satisfies_con3 <- con3(df$Var1, df$Var2) <= 0
df$satisfies_con4 <- con4(df$Var1, df$Var2) <= 0
df$satisfies_con5 <- con5(df$Var1, df$Var2) <= 0
df$feasible <- df$satisfies_con1 & df$satisfies_con2 & df$satisfies_con3 & df$satisfies_con4 & df$satisfies_con5

feasible_points <- df[df$feasible, c("Var1", "Var2")]

p<- plot_ly()
p <- add_trace(p, x = feasible_points$Var1, y = feasible_points$Var2, mode = "markers")
p


# Create constraint labels
annotations <- list(
  list(x = -4, y = con1(-4, 0), xref = "x", yref = "y", text = "x1 - 2x2 <= 0", showarrow = TRUE, arrowhead = 0, ax = -40, ay = -40),
  list(x = -4, y = con2(-4, 0), xref = "x", yref = "y", text = "2x1 - 3x2 <= 2", showarrow = TRUE, arrowhead = 0, ax = -40, ay = -40),
  list(x = -4, y = con3(-4, 0), xref = "x", yref = "y", text = "x1 - x2 <= 3", showarrow = TRUE, arrowhead = 0, ax = -40, ay = -40),
  list(x = -4, y = con4(-4, 0), xref = "x", yref = "y", text = "-x1 + 2x2 <= 2", showarrow = TRUE, arrowhead = 0, ax = -40, ay = -40),
  list(x = -4, y = con5(-4, 0), xref = "x", yref = "y", text = "-2x1 + x2 <= 0", showarrow = TRUE, arrowhead = 0, ax = -40, ay = -40)
)

# Define the shape of the feasible region
x_range <- range(df$Var1)
y_range <- range(df$Var2)
shapes <- list(
  list(
    type = "path",
    path = paste0("M", x_range[1], ",", con1(x_range[1], 0), "L", x_range[2], ",", con1(x_range[2],0), "L", x_range[2], ",", con4(x_range[2], y_range[2]), "L", x_range[1], ",", con5(x_range[1], y_range[2]), "Z"),
    fillcolor = "rgba(255, 0, 0, 0.2)",
    line = list(color = "transparent")
  )
)
                                                                                        
p <- layout(p, annotations = annotations, shapes = shapes)
p


    



theta0_true <- runif(1, 0, 1)
theta1_true <- runif(1, 0, 1)
N <- 1000
mean_error <- 0
sd_error <- 0.1
error <- rnorm(N, mean = mean_error, sd = sd_error)

x <- seq(from = 1, to = N, by = 1)/N
y <- theta0_true + theta1_true*x + error
regression <- lm(y ~ x)
# This fits a linear model to the data with y as the response variable and x as the predictor variable.

theta0_estimate <- coef(regression)[1]
theta1_estimate <- coef(regression)[2]
#This will give you the estimates of the true parameters based on the generated data, but again it is 
#important to keep in mind that this is a simulation and the estimates obtained may not reflect the true parameters.


estimate_parameters <- function(parameters, x, y) {
  theta0 <- parameters[1]
  theta1 <- parameters[2]
  y_predicted <- theta0 + theta1*x
  max_difference <- max(abs(y - y_predicted))
  return(max_difference)
}

initial_guess <- c(0, 0)
result <- optim(par = initial_guess, fn = estimate_parameters, x = x, y = y, method = "BFGS")
theta0_estimate <- result$par[1]
theta1_estimate <- result$par[2]



estimate_parameters <- function(parameters, x, y) {
  theta0 <- parameters[1]
  theta1 <- parameters[2]
  y_predicted <- theta0 + theta1*x
  sum_difference <- sum(abs(y - y_predicted))
  return(sum_difference)
}

initial_guess <- c(0, 0)
result <- optim(par = initial_guess, fn = estimate_parameters, x = x, y = y, method = "BFGS")
theta0_estimate <- result$par[1]
theta1_estimate <- result$par[2]

# Run 3 sets of experiments

#generate a sample of N = 1000 observation errors ξi from a uniform distribution on [-1,1]:
#uniform 
N <- 1000
error <- runif(N, min = -1, max = 1)
#gaussian
N <- 1000
error <- rnorm(N)
#cauchy
N <- 1000
error <- rcauchy(N)

summary(error)
