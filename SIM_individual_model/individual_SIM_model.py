import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

################################################## SIMULATION CONSTANTS

gamma = 1 / 14  # human recovery time
beta = 0.5  # transmission rate of malaria
mu = 1 / 14  # exposure period
sigma = 0  # immunity breakthrough rate

################################################## STATES

susceptible = 'S'
infected = 'I'
recovered = 'R'


################################################## INITIALISING GRAPH

def get_mosquito_transmission():

    # Uses the concept of vectorial capacity to determine the likelihood of a transmission between nodes
    a = np.random.uniform(low=0.1, high=2)       # rate at which a human is bitten by a mosquito
    b = np.random.uniform(low=0.5, high=0.5)     # proportion of infected bites that cause infection in the human host
    c = np.random.uniform(low=0.5, high=0.5)     # transmission efficiency from humans to mosquitoes
    m = np.random.uniform(low=0, high=10)        # number of mosquitoes in the region per human
    mu = np.random.uniform(low=0.14, high=0.23)  # life expectancy of mosquitoes

    # This will be used as a weight between nodes, giving the likelihood of transmission between humans along that node
    return (a ** 2 * b * c * m) / mu


def initialise_graph(num_nodes, num_infected_nodes):

    # Initialising the network structure
    G = nx.Graph()

    # Initialising the number of nodes
    all_nodes = np.arange(num_nodes)

    # Randomly choosing infected nodes
    random_infected_nodes = np.random.choice(all_nodes, size=num_infected_nodes, replace=False)

    # Looping through all nodes
    for node in all_nodes:

        # Checking if node has been randomly infected
        if np.any(node == random_infected_nodes):

            # Adding infected node if required
            G.add_node(node, state=infected, infects=[], infected_by=[])

        else:

            # Adding susceptible node if required
            G.add_node(node, state=susceptible, infects=[], infected_by=[])

    # Looping through edges
    for patch in G.nodes():
        for other_patch in G.nodes():

            # Connecting node to all other nodes except itself and setting distance between
            if patch != other_patch:
                G.add_edge(patch, other_patch, weight=get_mosquito_transmission())

    # Returning the resulting structure
    return G


################################################## RUNNING SIMULATION

def get_event_rates(G):
    # Stores the node that caused the event (e.g. I_k will cause S_j to become infected)
    source_nodes = []

    # Stores the node that the event happened to (S_j becoming infected due to I_k)
    target_nodes = []

    # Stores the propensity of the event
    event_rates = []

    # Looping through all nodes in the system
    for i, node in enumerate(G.nodes(data=True)):

        # Exctacting the node label and info (which includes the state, weight, etc)
        node_label, node_info = node

        if node_info['state'] == susceptible:

            # Determining the shortest path length (weight) to each nb of the current node
            all_nbs = nx.single_source_dijkstra_path_length(G, source=node_label, weight='weight')

            # Looping through all neighbours
            for nb in all_nbs:

                # Determining if the current neighbour is infected
                if G.nodes[nb]['state'] == infected:

                    # Determining the infected node's label and calculating its event rate based on its weight (mosquito dynamics)
                    infected_nb_label = nb
                    infected_nb_rate = all_nbs[nb]

                    # Storing the results
                    source_nodes.append(infected_nb_label)
                    target_nodes.append(node_label)
                    event_rates.append(infected_nb_rate)

        elif node_info['state'] == infected:

            # Storing the results
            source_nodes.append(node_label)
            target_nodes.append(node_label)
            event_rates.append(gamma)

    # Returning the results as numpy arrays
    return np.array(source_nodes), np.array(target_nodes), np.array(event_rates)


def get_chosen_state_position(event_rates, sum_rates):
    # Determining the probabilities of each event
    event_probabilities = event_rates / sum_rates

    # Randomly choosing an event from a discrete probability function
    random_position = np.random.choice(np.arange(len(event_probabilities)), size=1, p=event_probabilities)[0]

    # Returning the resulting position
    return random_position


def get_chosen_time(sum_rates):

    # Generating a uniform random number
    r1 = np.random.uniform()

    # Choosing current time by drawing from exponential
    chosen_time = np.log(1.0 / r1) / sum_rates

    # Returning the resulting time
    return chosen_time


def update_network(G, event_source, event_target):
    # Checking if the required state is susceptible
    if G.nodes()[event_target]['state'] == susceptible:

        # Setting state to infected
        G.nodes()[event_target]['state'] = infected

        # Storing infection information
        G.nodes()[event_source]['infects'].append(event_target)
        G.nodes()[event_target]['infects'].append(event_source)

    # Checking if the required state is infected
    elif G.nodes()[event_target]['state'] == infected:

        # Setting the state to recovered
        G.nodes()[event_target]['state'] = recovered

    # Returning the network system
    return G


def run_gillespie(G, S0, I0, R0, t_max):

    # Initialising time
    t = 0

    # Initialising results to store compartment totals
    S_data, I_data, R_data, t_data = [S0], [I0], [R0], [t]

    # Looping through time until convergence or maximum simulation time is reached
    while t < t_max:

        # Determining the possible events and their rates
        source_nodes, target_nodes, event_rates = get_event_rates(G)
        sum_rates = np.sum(event_rates)

        # Checking for convergence
        if len(event_rates) == 0:
            break

        # Choosing the event to occur by drawing from a discrete probability distribution
        chosen_position = get_chosen_state_position(event_rates, sum_rates)

        # Determining which event took place
        event_source = source_nodes[chosen_position]  # node causing the event (e.g. I_j causes an infection to S_k)
        event_target = target_nodes[chosen_position]  # node event is occurring to (e.g. S_k becomes infected by I_j)

        # Determining the time for the event to occur
        chosen_time = get_chosen_time(sum_rates)

        # Updating the time
        t += chosen_time

        # Updating the network
        G = update_network(G, event_source, event_target)

        # Saving the results
        S_data.append(np.sum([G.nodes()[node]['state'] == 'S' for node in G.nodes()]))
        I_data.append(np.sum([G.nodes()[node]['state'] == 'I' for node in G.nodes()]))
        R_data.append(np.sum([G.nodes()[node]['state'] == 'R' for node in G.nodes()]))
        t_data.append(t)

    return G, S_data, I_data, R_data, t_data


################################################## DISPLAYING GRAPH

def show_graph(G):

    # Setting colours for the different node types
    color_map = []

    # Looping through all possible nodes and extracting their label and data (which includes their state)
    for node_label, node_data in G.nodes(data=True):

        # Setting susceptible nodes to blue
        if node_data['state'] == susceptible:
            color_map.append('tab:blue')

        # Setting infected nodes to red
        elif node_data['state'] == infected:
            color_map.append('tab:red')

        # Setting immune nodes to green
        else:
            color_map.append('tab:green')

    # Drawing the graph
    plt.clf()
    plt.figure()
    pos = nx.spring_layout(G)
    nx.draw_networkx(G, pos, node_size=100, node_color=color_map, with_labels=False)
    # plt.legend()
    plt.axis('off')
    plt.show()


################################################## CREATING AN INFECTION TREE

def get_infection_tree(G):

    # Creating infection tree as a directed graph
    infection_tree = nx.DiGraph()

    # Adding nodes
    infection_tree.add_nodes_from(range(G.number_of_nodes()))

    # Adding edges (i.e. who infected who)
    for node in G.nodes(data=True):
        node_label, node_info = node

        # Checking if the current node has caused any infections
        if len(node_info['infects']) > 1:

            # Adding the nodes infected by the current node to the network
            for infected in node_info['infects']:
                infection_tree.add_edge(node_label, infected)

    # Removing isolated nodes
    infection_tree.remove_nodes_from(list(nx.isolates(infection_tree)))

    # Returning the tree structure
    return infection_tree


################################################## MAIN PROGRAM

# Network parameters
num_nodes = 100

# Initial conditions
I0, R0 = 1, 0
S0 = num_nodes - I0 - R0

# Initialising the network
G = initialise_graph(num_nodes, I0)

# Setting maximum simulation time
max_time = 100

# Running the algorithm
G, S_data, I_data, R_data, t_data = run_gillespie(G, S0, I0, R0, max_time)
infection_tree = get_infection_tree(G)

# Plotting the results
plt.clf()
plt.figure()
plt.plot(t_data, S_data, label='Susceptible')
plt.plot(t_data, I_data, label='Infected')
plt.plot(t_data, R_data, label='Immune')
plt.title('Evolution of compartment sizes over time for an individual-based model \nwith one variant')
plt.xlabel('Time (days)')
plt.ylabel('Compartment Size')
plt.legend()
plt.savefig('individual_SIM_CompartmentalData.png')

# Drawing the family tree
plt.clf()
plt.figure()
pos = nx.spring_layout(G)
nx.draw(infection_tree, pos, with_labels=True, arrowsize=10)
plt.savefig('individual_SIM_InfectionTree.png')
