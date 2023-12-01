# MPhys_Project-Epidemic_Modelling
## Compartmental Models
Both SIM (susceptible-infected-immune) and SEIM (susceptible-exposed-infected-immune) compartmental models were used to simulate the impact of catastrophic events on the transmission dynamics of a vector-borne disease. The models consisted of ordinary differential equations, and used a time-dependent death rate to incorporate shock events with the intent of producing population bottlenecks. The mosquito-based transmission dynamics were encapsulated into the probability of transmission between the compartments, as opposed to having an explicit mosquito subset, allowing for a simple model to be analysed. The set of ordinary differential equations were numerically-approximated using the Euler method.

## Individual-Based Model
An individual-based model was implemented to simulate the transmission dynamics of a vector-borne disease within a susceptible, infected, and immune population. The model is a fully-connected network, where each node represents a human individual. The edges between the nodes represent the likelihood of an infection being spread by mosquitoes between the two nodes. This likelihood was calculated by considering the capacity of the vectors to transmit disease, according to the formula
$$p_{transmission} = \frac{a^2bcm}{\mu}$$
where $a$ is the rate at which a single human is bitten by the mosquitoes, $b$ is the proportion of bites from an infected mosquito that cause an infection to develop within the human, $c$ is the transmission efficiency of the disease from humans to mosquitoes, $m$ is the ratio of mosquitoes to humans at the location, and $\mu$ is the life expectancy of a mosquito.

The Gillespie algorithm was used to simulate the dynamics within the network, allowing for increased computational efficiency and a stochastic simulation. The simulation was run until convergence - meaning all individuals recovered - or until the maximum simulation time was reached. During each iteration of the simulation, the following algorithm was completed:
* Initialise the network with the initial conditions for $S_0$, $I_0$, and $M_0$.
* Determine each possible event $X_i$. The possible events include an infected node recovering, or a susceptible node becoming infected. Therefore, for a given susceptible node, all infected neighbours - and the transmission probabilities along the edges - were determined.
* Calculate each event rate $a_i$. The rates may change and so must be recalculated each iteration.
* Calculate the sum of the rates $a = \sum_{i=1}a_i$.
* Determine the duration of the event by selecting a time step $dt$ from $dt = \frac{1}{a}log(\frac{1}{r1})$, where $r1$ is drawn from a uniform distribution.
* Determine which event will occur with probability $P(X_i)=\frac{a_i}{a}$.
* Perform the chosen event and update the network.
* Update the simulation time by $t \to t + dt$.
* Check for convergence (meaning no more infected nodes are present in system) or if the maximum simulation time has been reached. If so, stop the simulation, and otherwise return to Step 2.
