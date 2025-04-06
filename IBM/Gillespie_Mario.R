################################################################################
#
# Demographic Functions
# 
#
#
################################################################################



# Required Packages
library("sads")

# Define Demographic Functions

# Birth: Givin the abundances of taxa, randomly pick one taxa and add one unit
birth <- function(N = N, b = b, otus = otus){
  p_i = b*N
  i = sample(seq_along(otus), size = 1, prob = p_i)
  N[i] <<- N[i] + 1
}

# Death: Givin the abundances of taxa, randomly pick one and remove one unit
#         The unit removed will be added to a "relic" pool
death <- function(N = N, R = R, d = d, otus = otus){
  p_i = d*N
  i = sample(seq_along(otus), size = 1, prob = p_i)
  N[i] <<- N[i] - 1
  R[i] <<- R[i] + 1
}

# Decay: Units in the relic pool degarde over time. 
#        There are a few ways this could happen
#         1. Random decay - uniform decay based on abundance
#         2. Species specific decay - some species decay faster and some slower  
#            which we model using a beta distribution 
#         3. Density dependent decay - species with more individuals in the 
#            relic pool will decay faster. Why, well more abundent taxa are 
#            close together and therefore may degrade their own relic DNA. 
decay <- function(R = R, u = u, otus = otus, 
                  a.decay = 0.7, b.decay = 0.7,
                  method = c("random", "sp.spec", "den.depend")){
  if (method == "random"){
    sp.decay <- rep(1, length(otus))
    d.prob <- (sp.decay * (R/sum(R)))
  } 
  if (method == "sp.spec"){
    sp.decay <- round(rbeta(length(otus), shape1 = a.decay, shape2 = b.decay), 3)
    # hist(rbeta(length(otus), shape1 = a.decay, shape2 = b.decay))
    d.prob <- (sp.decay * (R/sum(R)))
    } 
  if (method == "den.depend"){
    sp.decay <- rep(1, length(otus))
    d.prob <- ((sp.decay * (R/sum(R)))^2) + 1/200 * (sp.decay * (R/sum(R)))
    }
  if (method %in% c("random", "sp.spec", "den.depend") == FALSE){
    stop("You must select one of the following methods: random, sp.spec, or den.depend")
  }
  p_i = u*d.prob
  i = sample(seq_along(otus), size = 1, prob = p_i)
  R[i] <<- R[i] - 1
}

# Immigration: to prevnt taxa from going extinct, we need new individuals coming 
#              in from the regional species pool 
immigration <- function(N = N, G = G, im = im, otus = otus){
  p_i = im*G
  i = sample(seq_along(otus), size = 1, prob = p_i)
  N[i] <<- N[i] + 1 
}

# Emmigration: it is also possible for individuals to leave the local community
emigration <- function(N = N, em = em, otus = otus){
  p_i = em*N
  i = sample(seq_along(otus), size = 1, prob = p_i)
  N[i] <<- N[i] - 1 
}

# Zombie: In addition to the immigration of intact species, we can also have the 
#         immigration of relic DNA. We call this zombie immigration
zombie <- function(R = R, G = G, z = z, otus = otus){
  p_i = z*G
  i = sample(seq_along(otus), size = 1, prob = p_i)
  R[i] <<- R[i] + 1 
}

# Define Simulation Model
IR_model <- function(steps = "", b = "", d = "", u = "", z = "",
                     im.b = "", em.b = "", init = "",
                     decay.m = c("random", "sp.spec", "den.depend"), 
                     a.decay = "", b.decay = ""){
 
  # Define Initial Simulation Options
  num.steps <- steps
  sim.time <- 0 # reset initial time
  n.birth <- 0; n.death <- 0; n.decay <- 0; 
  n.immigration <- 0; n.emigration <- 0; n.zombie <- 0

  # Define OTUs and Rel Abundance in Regional Pool
  otus <- paste("OTU", sprintf("%05d", seq(1:200)), sep = "")
  #G <- rlnorm(n=length(otus), meanlog = 1, sdlog = 0.98)
  G <- rls(n = 1:length(otus), N = length(otus), alpha = 5)
  G <<- G[rev(order(G))] # / sum(G)
  #plot(G)
  #plot(log10(G))

  # Initialize Local Community and Relic Tables 
  N <<- table(otus) * 0; n <- 0
  R <<- table(otus) * 0; r <- 0

  # Initialize Local Community Membership
  comm <- sample(otus, size = 100, replace = T, prob = G)
  ind <- which(otus %in% comm)
  N[ind] <<- N[ind] + table(comm)

  # Save Time and Events 
  ints <- vector(mode = "numeric", length = 0)
  event.list <- vector(mode = "character", length = 0)
  N_series <- vector(mode = "numeric", length = 0)
  R_series <- vector(mode = "numeric", length = 0)
  T_series <- vector(mode = "numeric", length = 0)
  
  # Run Simulation
  while(sim.time < num.steps){
    
    n <- sum(N)
    r <- sum(R)

    im <- if (n < init) {((init + im.b) - n)} else {im.b}; 
    em <- if (n < init) {floor(0.1 * n)} else {em.b}; 
  
    # Timestep
    interval <- rexp(1, rate = n*b + n*d + r*u + im + em + z)
    ints <- c(ints, interval)
    sim.time <- sim.time + interval
  
    # Pick Event
    event = sample(c('birth', 'death', 'decay', 
                     'immigration', 'emigration', 'zombie'),
                   size = 1, prob = c(b*n, d*n, r*u, im, em, z))
    
    switch(event,
           birth = birth(N = N, b = b, otus = otus),
           death = death(N = N, R = R, d = d, otus = otus),
           decay = decay(R = R, u = u, otus = otus, 
                         method = decay.m, 
                         a.decay = a.decay, b.decay = b.decay),
           immigration = immigration(N = N, G = G, im = im, otus = otus),
           emigration = emigration(N = N, em = em, otus = otus),
           zombie = zombie(R = R, G = G, z = z, otus = otus),
           stop("event must be either \"birth\", \"death\", \"decay\",
                \"immigration\", \"emigration\", or \"zombie\""))
    
    # Store Event
    event.list <- c(event.list, event)
    N_series <- c(N_series, sum(N))
    R_series <- c(R_series, sum(R))
    T_series <- c(T_series, sim.time)
    
    # Print Statement
    print(paste("Time = ", round(sim.time, 3), ": ", event, 
                "; N = ", sum(N), "; R = ", sum(R), sep = ""))
  }
  
  # Summary and Print Statements
  event_summary <-  table(Events = event.list)
  time_series <- data.frame(Time = T_series, N = N_series, R = R_series)
  print(event_summary)
  print(data.frame(probs = cbind(birth = b*n, death = d*n, decay = r*u, 
                       immigration = im, emigration = em, zombie = z)))
  print(data.frame(state = cbind(birth = b, death = d, decay = u, 
                      immigrtion = im, emigration = em, zombie = z)))
  return(list(N = N, R = R, intervals = ints, 
              events = event.list, time_series = time_series))
}

# IR_sim <- function()

# Equal Birth and Death and Decay
model1 <- IR_model(steps = 20, b = 0.2, d = 0.2, u = 0.2, z = 10,
                   im.b = 50, em.b = 10, init = 500,
                   decay.m = "random")

# Rank Abundance Plot
K <- model1$N + model1$R
plot(rad(as.numeric(K) / sum(K)), type = "l", col = rgb(0,0,0,0.5), lwd = 5)
lines(rad(as.numeric(N) / sum(N)), col = rgb(0,0,1,0.5), lwd = 5)
lines(rad(as.numeric(R) / sum(R)), col = rgb(1,0,0,0.5), lwd = 5)
legend("topright", legend = c("Total", "Intact", "Relic"), bty = "n",
       pch = 22, col = c("black", "blue", "red"), 
       pt.bg = c(rgb(0,0,0,0.5), rgb(0,0,1,0.5), rgb(1,0,0,0.5)))

ks.test(model1$N, model1$R)

# Cumulative Abundance Distribution
# nN <- seq(1, max(N), 1)
# cN <- rep(NA, length(nN))
# for(i in 1:length(cN)){
#   cN[i] <- sum(N[which(N > i)])
# }
# cNrel <- cN / sum(N)
# 
# nR <- seq(1, max(R), 1)
# cR <- rep(NA, length(nR))
# for(i in 1:length(cR)){
#   cR[i] <- sum(R[which(R <= i)])
# }
# cRrel <- cR / sum(R)
# 
# plot(log10(cNrel) ~ log10(nN), col = "blue")
# points(cRrel ~ nR, col = "red")
# max(R)
# tR <- table(R)

replicate(10, {
  model1 <- IR_model(steps = 20, b = 0.2, d = 0.2, u = 0.2, z = 10,
                     im.b = 50, em.b = 10, init = 500,
                     decay.m = "random")
  
  # Rank Abundance Plot
  K <- model1$N + model1$R
  lines(rad(as.numeric(K) / sum(K)), col = rgb(0,0,0,0.5), lwd = 5)
  lines(rad(as.numeric(N) / sum(N)), col = rgb(0,0,1,0.5), lwd = 5)
  lines(rad(as.numeric(R) / sum(R)), col = rgb(1,0,0,0.5), lwd = 5)
  ks.test(model1$N, model1$R)
})

# Slow Decay
rm(list = c("G", "N", "R"))
model2 <- IR_model(steps = 200, b = 0.2, d = 0.2, u = 0.02, z = 10,
                   im.b = 50, em.b = 10, init = 500,
                   decay.m = "random")

# Rank Abundance Plot
plot(rad(as.numeric(G)), type = "l")
K <- N + R
plot(rad(as.numeric(K) / sum(K)), type = "l", col = rgb(0,0,0,0.5), lwd = 5)
lines(rad(as.numeric(N) / sum(N)), col = rgb(0,0,1,0.5), lwd = 5)
lines(rad(as.numeric(R) / sum(R)), col = rgb(1,0,0,0.5), lwd = 5)

ks.test(model2$N, model2$R)

# Very Slow Decay
rm(list = c("G", "N", "R"))
model2 <- IR_model(steps = 500, b = 0.2, d = 0.2, u = 0.0002, z = 10,
                   im.b = 50, em.b = 10, init = 500,
                   decay.m = "random")

# Rank Abundance Plot
K <- N + R
plot(rad(as.numeric(K) / sum(K)), type = "l", col = rgb(0,0,0,0.5), lwd = 5)
lines(rad(as.numeric(N) / sum(N)), col = rgb(0,0,1,0.5), lwd = 5)
lines(rad(as.numeric(R) / sum(R)), col = rgb(1,0,0,0.5), lwd = 5)

replicate(10, {
  model2 <- IR_model(steps = 500, b = 0.2, d = 0.2, u = 0.0002, z = 10,
                     im.b = 50, em.b = 10, init = 500,
                     decay.m = c("random", "sp.spec", "den.depend"), 
                     a.decay = 0.7, b.decay = 0.7)
  
  # Rank Abundance Plot
  K <- N + R
  lines(rad(as.numeric(K) / sum(K)), col = rgb(0,0,0,0.5), lwd = 5)
  lines(rad(as.numeric(N) / sum(N)), col = rgb(0,0,1,0.5), lwd = 5)
  lines(rad(as.numeric(R) / sum(R)), col = rgb(1,0,0,0.5), lwd = 5)
  
  ks.test(model2$N, model2$R)
})



sum(N); sum(R)
K <- N + R; sum(K)
sum(N) / sum(R)

# Sample then test
N.samp <- table(sample(names(model2$N), size = 1000, replace = T, prob = model2$N))
R.samp <- table(sample(names(model2$R), size = 1000, replace = T, prob = model2$R))

ks.test(N.samp, R.samp)





# Initialize Local Community Membership
comm <- sample(otus, size = 100, replace = T, prob = G)
ind <- which(otus %in% comm)
N[ind] <<- N[ind] + table(comm)

plot(0, 0, xlim = c(0, 200), ylim = log10(c(1, 25000)), type = "n",
     xlab = "Time", ylab = "Populaiton Size (log10)")
legend("topleft", c("Intact", "Relic"), pch = 16, col = c("blue", "red"), bty = "n")

points(x = sim.time, y = log10(sum(N)), pch = 16, cex = 0.2, col = "blue")
points(x = sim.time, y = log10(sum(R)), pch = 16, cex = 0.2, col = "red")


  
