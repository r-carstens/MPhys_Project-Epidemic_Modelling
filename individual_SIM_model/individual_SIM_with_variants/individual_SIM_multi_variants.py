import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd


################################################## SIMULATION CONSTANTS

gammas = [1/14, 1/10] # human recovery times for variants
sigmas = [0, 0]       # immunity breakthrough rates from variants

connection_prob = 0.3 # probability that a node is connected to another node


################################################## STATES

susceptible = 'S'
infected = 'I'
recovered = 'R'


################################################## INITIALISING GRAPH

def get_mosquito_transmission_weights():

    # Uses the concept of vectorial capacity to determine the likelihood of a transmission between nodes
    a = np.random.uniform(low=0.1, high=2)       # rate at which a human is bitten by a mosquito
    b = np.random.uniform(low=0.5, high=0.5)     # proportion of infected bites that cause infection in the human host
    c = np.random.uniform(low=0.5, high=0.5)     # transmission efficiency from humans to mosquitoes
    m = np.random.uniform(low=0, high=5)        # number of mosquitoes in the region per human
    mu = np.random.uniform(low=0.14, high=0.23)  # life expectancy of mosquitoes

    # This will be used as a weight between nodes, giving the likelihood of transmission between humans along that node
    return (a ** 2 * b * c * m) / mu


def get_graph_nodes(G, all_nodes, num_variants):
    
    # Randomly choosing infected nodes (currently only 1 infection per variant in the system)
    random_infected_nodes = np.random.choice(all_nodes, size=num_variants, replace=True)
    
    for node in all_nodes:
        
        # Creating a structure to store each node's states with respect to the variants
        initial_states = np.array([susceptible for i in range(num_variants)])
        
        # Checking if the node is infected with any variant and altering the initial states to include the infection(s)
        if len(np.where(node == random_infected_nodes)[0]) > 0:
            initial_states[np.where(node == random_infected_nodes)] = infected
        
        # Adding the node
        G.add_node(node, states=initial_states)
        
    return G


def get_graph_edges(G):
    
    # Looping through all possible edges
    for patch in G.nodes():
        for other_patch in G.nodes():

            # Connecting node to all other nodes except itself and setting distance between
            if patch != other_patch and np.random.uniform() < connection_prob:
                G.add_edge(patch, other_patch, weight=get_mosquito_transmission_weights())
                
    return G


def initialise_graph(num_nodes, num_variants):

    # Initialising the network structure and the node labels
    G = nx.Graph()
    all_nodes = np.arange(num_nodes)
    
    # Creating the graph (these are functions to allow for more graphs to be easily constructured)
    G = get_graph_nodes(G, all_nodes, num_variants)
    G = get_graph_edges(G)
                    
    return G


################################################## RUNNING SIMULATION

def get_variant_populations(G, variant):
    
    # Determining all nodes susceptible to and infected by the current variant
    susceptible_nodes = [node_label for node_label, node_data in G.nodes(data=True) if susceptible in node_data['states'][variant]]
    infected_nodes = [node_label for node_label, node_data in G.nodes(data=True) if infected in node_data['states'][variant]]

    return susceptible_nodes, infected_nodes
    
    
def get_variant_recovery_data(G, infected_nodes, variant):
    
    # Determining potential recoveries
    recovery_rates = [gammas[variant] for infected_node in infected_nodes]
    
    return recovery_rates


def get_variant_infection_data(G, susceptible_nodes, infected_nodes, variant):
    
    sources, targets, infection_rates = [], [], []
    
    # Determining potential infections
    for sus_node in susceptible_nodes:
            
        # Finding infected neighbours and their transmission probabilities
        current_infected_nbs = set(G.neighbors(sus_node)) & set(infected_nodes)
        current_infection_rates = [G.get_edge_data(sus_node, infected_nb)['weight'] for infected_nb in current_infected_nbs]
        
        # Storing the results
        sources += current_infected_nbs
        targets += [sus_node for i in range(len(current_infected_nbs))]
        infection_rates += current_infection_rates
        
    return sources, targets, infection_rates


def get_event_rates(G, num_variants):
    
    # Data to store event results
    sources, targets, rates, variants = [], [], [], []
    
    # Determining event data for all variants
    for variant in range(num_variants):

        # Subgrouping into susceptible and infected nodes for the current variant
        susceptible_nodes, infected_nodes = get_variant_populations(G, variant)
        
        # Determining all possible event data for the current variant
        variant_recovery_rates = get_variant_recovery_data(G, infected_nodes, variant)
        variant_infected_sources, variant_susceptible_targets, variant_infection_rates = get_variant_infection_data(G, susceptible_nodes, infected_nodes, variant)

        # Storing the results
        sources += (infected_nodes + variant_infected_sources)
        targets += (infected_nodes + variant_susceptible_targets)
        rates += (variant_recovery_rates + variant_infection_rates)
        variants += [variant for i in range(len(variant_recovery_rates) + len(variant_infection_rates))]
        
    return sources, targets, rates, variants
    
    
def get_chosen_state_position(event_rates, sum_rates):
    
    # Determining the probabilities of each event and randomly choosing an event from a discrete probability function
    event_probabilities = event_rates / sum_rates    
    random_position = np.random.choice(np.arange(len(event_probabilities)), size=1, p=event_probabilities)[0]
    
    return random_position


def get_chosen_time(sum_rates):
    
    # Generating a uniform random number and choosing timestep by drawing from exponential
    r1 = np.random.uniform()
    chosen_time = np.log(1.0 / r1) / sum_rates
    
    return chosen_time


def update_network(G, event_source, event_target, event_variant):
            
    # Checking if the required variant state is susceptible
    if G.nodes()[event_target]['states'][event_variant] == susceptible:
        
        # Setting to infected
        G.nodes()[event_target]['states'][event_variant] = infected
        
    # Checking if the required variant state is infected
    elif G.nodes()[event_target]['states'][event_variant] == infected:
        
        # Setting the state to recovered
        G.nodes()[event_target]['states'][event_variant] = recovered
                
    return G


def get_current_totals(G, num_variants):
    
    # Creating structures to store the current totals
    current_data = np.zeros(shape=(num_variants, 3))
    
    # Looping through each possible variant
    for variant in range(num_variants):
        
        # Determining the population totals
        num_susceptible = np.sum([G.nodes()[node]['states'][variant] == susceptible for node in G.nodes()])
        num_infected = np.sum([G.nodes()[node]['states'][variant] == infected for node in G.nodes()])
        num_recovered = np.sum([G.nodes()[node]['states'][variant] == recovered for node in G.nodes()])
        
        # Storing the current results
        current_data[variant] = num_susceptible, num_infected, num_recovered
    
    return current_data.reshape(num_variants * 3)
                                 
                                   
def run_gillespie(G, num_variants, t_max):
    
    # Creating file to store data
    out_file = open('individual_SIM_variants_data.txt', 'w')
    out_file.write('Time,' + ','.join([f'S_{variant},I_{variant},R_{variant}' for variant in range(num_variants)]) + '\n')
    
    # Setting initial conditions
    t = 0
    
    while t < t_max:
    
        # Determining the possible events and their rates
        sources, targets, rates, variants = get_event_rates(G, num_variants)
        sum_rates = np.sum(rates)
        
        # Checking for convergence
        if len(rates) == 0:
            break

        # Choosing the event to occur by drawing from a discrete probability distribution
        chosen_position = get_chosen_state_position(rates, sum_rates)

        # Determining which event took place
        event_source = sources[chosen_position]    # node causing the event (e.g. I_j causes an infection to S_k)
        event_target = targets[chosen_position]    # node event is occurring to (e.g. S_k becomes infected)
        event_variant = variants[chosen_position]  # specific variant involved
        
        # Updating the system
        G = update_network(G, event_source, event_target, event_variant)
        t += get_chosen_time(sum_rates)

        # Storing the results
        population_results = get_current_totals(G, num_variants)
        out_file.write(str(t) + ',' + ','.join(population_results.astype(str)) + '\n')

    out_file.close()
    return G


# Initialising the system
N = 1000
num_variants = 2

# Setting the initial conditions
I0, R0 = 1, 0
S0 = N - I0 - R0
t_max = 100

# Creating the graph and running the simulation
G = initialise_graph(N, num_variants)
G = run_gillespie(G, num_variants, t_max)

# Reading the data
data = pd.read_csv('individual_SIM_variants_data.txt')
data.columns

# Plotting the data
for i in range(1, len(data.columns)):
    plt.plot(data[data.columns[0]], data[data.columns[i]], label=data.columns[i])
    
plt.legend()
plt.title('Evolution of population sizes over time for an individual-based model \nwith %s variants' % num_variants)
plt.xlabel('Time (days)')
plt.ylabel('Population Size')
plt.legend()
plt.savefig('individual_SIM_variants_results.png')
