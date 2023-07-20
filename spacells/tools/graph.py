import networkx as nx
import matplotlib.pyplot as plt



class graph:
    """[summary]
        create graph for the target cluster from input ndarray of points
    """
    def __init__(self, ndarray_in):
        self.processed_input  = []
        self.G = nx.Graph()
        self.position={}
        self.create_graph(ndarray_in)

    def ndarray_to_tuple(self, ndarray_in):
        """[summary]

        Args:
            ndarray_in (numpy.array): input array for points
        """
        self.processed_input = [[tuple(ie) for ie in e] for e in ndarray_in]


    def create_graph(self, ndarray_in):
        """[summary]

        Args:
            ndarray_in (numpy.array): input array for points
        """
        self.ndarray_to_tuple(ndarray_in)
        for i in self.processed_input:
            nx.add_cycle(self.G, i)
            for point in i:
                self.position[point] = point

    def get_graph_info(self):
        """[summary] 
            retrieve information about the created graph
        """
        print("G nodes: ", self.G.nodes()) 
        print("G edges: ", self.G.edges()) 
        print("G diameter: ", nx.diameter(self.G))

    def get_graph_diameter(self):
        """[summary]

        Returns:
            int: diameter of the created graph with unit length, 1 for each edge
        """
        return nx.diameter(self.G)    

    def get_graph(self):
        """[summary]

        Returns:
            Graph: the graph of the cluster
        """
        return self.G

    def draw_graph(self, output=None):
        """[summary]
            draw and save the graph in jpg format
        """

        fig, ax = plt.subplots()
        nx.draw(self.G, pos=self.position, node_color='k', ax=ax)
        nx.draw(self.G, pos=self.position, node_size=100, ax=ax)  # draw nodes and edges
        nx.draw_networkx_labels(self.G, pos=self.position)  # draw node labels/names

        # draw edge weights
        labels = nx.get_edge_attributes(self.G, 'weight')
        nx.draw_networkx_edge_labels(self.G, self.position, edge_labels=labels, ax=ax)
        
        plt.axis("on")
        xmin, xmax = ax.get_xlim()
        ax.set_xlim(xmin-0.25, xmax+0.25)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin-0.25, ymax+0.25)
        ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)

        if output != None:
            plt.savefig(output)  #'cluster.jpg'
