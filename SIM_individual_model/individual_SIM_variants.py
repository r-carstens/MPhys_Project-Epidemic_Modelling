import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


################################################## SIMULATION CONSTANTS

gamma = 1/14          # human recovery time from variant 1
sigma = 0             # immunity breakthrough rate from variant 1

gamma = 1/10          # human recovery time from variant 2
sigma = 0             # immunity breakthrough rate from variant 2

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
    m = np.random.uniform(low=0, high=10)        # number of mosquitoes in the region per human
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
                    
    # Returning the resulting graph
    return G

    
num_nodes = 100
num_variants = 2

G = initialise_graph(num_nodes, num_variants)
get_event_rates(G)
