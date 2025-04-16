import numpy as np
import C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset as simplex2
from distinct_permutations import distinct_permutations
from itertools import chain, combinations, product
from tools import sort_np_array_rows_lexicographically

#returns powerset without the initial empty tupel
def powerset(s: list):
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))

def recursive_step(l: list):
    new_l = l.copy()
    original_l = len(new_l)
    for i in range(original_l):
        if new_l[i][-1].count(0) == 1:
            return new_l
        else:
            for j in range(len(new_l[i][-1])):
                if(new_l[i][-1][j] == 0):
                    copy = new_l[i][-1].copy()
                    copy[j] += 1
                    edit = new_l[i].copy()
                    edit.append(copy)
                    new_l.append(edit)
    for i in range(original_l):
        new_l.pop(0)
    return recursive_step(new_l)

class Vertex(simplex2.Eji_LinComb):

    #Note that adding does not commute with other class operations, in particular not with value or to_array
    #e.g. a.value(dim)+b.value(dim) != (a+b).value(dim)
    def __add__(self, other):
        ans = Vertex.__new__(Vertex)
        ans._index = self._index+other._index
        ans._eji_counts = self._eji_counts+other._eji_counts
        return ans

    def value(self, dim):
        return self.get_canonical_form().hash_to_point_in_unit_hypercube(dim)

    def to_array(self):
        return self.get_canonical_form()._eji_counts

    def get_canonical_form(self):
        ans = Vertex.__new__(Vertex)
        ans._index = self._index
        ans._eji_counts = sort_np_array_rows_lexicographically(self._eji_counts)
        return ans
        
class Simplex():
    
    def __init__(self, n: int, k: int, vlist: list[Vertex]):
        self.n = n
        self.k = k
        self.dim = 2*(n-1)*k+1
        self.vlist = vlist
        self.num = len(vlist)

    def __eq__(self, other):
        if (self.n != other.n) or (self.k!=other.k) or (self.num!=other.num): return False
        array_list_self = []
        array_list_other = []
        for i in range(self.num):
            array_list_self.append(self.vlist[i].to_array())
            array_list_other.append(other.vlist[i].to_array())
        a = sorted(array_list_self, key = lambda x: tuple(x.flatten()))
        b = sorted(array_list_other, key = lambda x: tuple(x.flatten()))
        return all(np.array_equal(x, y) for x,y in zip(a, b))

    def __add__(self, other):
        assert(self.n == other.n and self.k == other.k)
        temp = self.vlist.copy()
        return Simplex(self.n , self.k, self.vlist+other.vlist)

    def get_canonical_form(self):
        new_vertices = [vertex.get_canonical_form() for vertex in self.vlist]
        return Simplex(self.n, self.k, new_vertices)
                

    def barycentric_subdivision(self):
        barycentre = sum(self.vlist, start=Vertex(self.n, self.k))
        simplex_list =[]
        seed_list = []
        a = np.zeros(self.num)
        a[0] += 1
        for perm in distinct_permutations(a):
            seed_list.insert(0, [list(perm)])
        allowed_vertex_combos = recursive_step(seed_list)
        for combo in allowed_vertex_combos:
            temp_vlist = [barycentre]
            for ele in combo:
                new_vertex = sum([self.vlist[i] for i in filter(lambda i: ele[i], range(self.num))], start = Vertex(self.n,self.k))
                temp_vlist.append(new_vertex)
            simplex_list.append(Simplex(self.n, self.k, temp_vlist))
        return simplex_list
    

        
    def projected_point(self, p: np.ndarray):
        assert len(p) == self.dim
        n = self.num
        b = np.zeros(n)
        A = np.zeros(shape=(n,n))
        for i in range(n):
            b[i] = np.dot(p, self.vlist[i].value(self.dim))
            for j in range(n):
                A[i,j] = np.dot(self.vlist[i].value(self.dim), self.vlist[j].value(self.dim))
        lamda  = np.linalg.solve(A, b)
        projection = np.zeros(self.dim)
        for i in range(n):
            projection += lamda[i]*self.vlist[i].value(self.dim)
        return projection, lamda, self.vlist

    def distance_to_point(self, p: np.ndarray):
        projection = self.projected_point(p)[0]
        return np.linalg.norm(p-projection, ord=2)


class SimplexMap():
    
    def __init__(self, n, k, subdivided = True):
        self.n = n
        self.k = k
        temp_list = self.simplices_across_all_k(self.n, self.k)
        print("Pre_subdivision length: ", len(temp_list))
        if subdivided:
            simplex_list = []
            for simplex in temp_list:
                simplex_list += simplex.barycentric_subdivision()
        else:
            simplex_list = temp_list
        print("Post_subdivision length: ", len(simplex_list))
        self.slist = self.remove_equivalent_simplices(n, k, simplex_list)
        print([s.vlist for s in self.slist])
        print("Post_reduction length: ", len(self.slist))


    def list_to_vertex(self, lt: list, n, k, kmax):
        assert len(lt) == n
        vertex = Vertex(n,kmax,[])
        for l in range (n):
            if(lt[l] != 0):
                vertex.add(simplex2.Maximal_Simplex_Vertex({simplex2.Eji(j=l,i=k)}))
        return vertex

    def generate_simplices_for_single_k(self, n, k, kmax):
        seed_list = []
        a = np.zeros(n)
        a[0] += 1
        for perm in distinct_permutations(a):
            seed_list.insert(0, [list(perm)])
        temp_list = recursive_step(seed_list)
        simplex_list = []
        for ele in temp_list:
            vlist = []
            for v in ele:
                vlist.append(self.list_to_vertex(v, n, k, kmax))
            simplex_list.append(Simplex(n, kmax, vlist))
        return simplex_list


    def simplices_across_all_k(self, n, k):
        superlist = []
        simplex_list = []
        for i in range(k):
            superlist.append(self.generate_simplices_for_single_k(n, i, k))
        for simplex_comb in product(*superlist):
            simplex_list.append(sum(simplex_comb, start=Simplex(n,k, [])))
        return simplex_list

    def remove_equivalent_simplices(self, n ,k, simplex_list):
        reduced_list = []
        for simplex in simplex_list:
            if simplex not in reduced_list:
                reduced_list.append(simplex.get_canonical_form())
        return reduced_list
        
    
    def choose_simplex(self, p: np.ndarray):
        dist_list = []
        simplex_list = self.slist
        print("Default: ", simplex_list[0].vlist)
        dist_list = [x.distance_to_point(p) for x in simplex_list]
        chosen_simplex = simplex_list[min(enumerate(dist_list), key = lambda x: x[1])[0]]
        print("Chosen: ", chosen_simplex.vlist)
        return chosen_simplex

    def get_lin_comb(self, p: np.ndarray):
        simplex = self.choose_simplex(p)
        lin_comb = simplex.projected_point(p)[1]
        return lin_comb, simplex.vlist


class Decoder():                

    def decode(self, n ,k , data):
        ans = np.zeros(shape = (n,k))
        for i in range(n):
            ans[i] += data[:k]

        data = data[k:]
        
        deltas, ejis = SimplexMap(n,k).get_lin_comb(data)
        ejis = [eji.to_array() for eji in ejis]
        print(deltas)
        for i in range(len(deltas)):
            ans += deltas[i]*ejis[i]
        
        return ans

if __name__ == "__main__":
    n = 2
    k = 2
    some_input = np.asarray([[4, 2],[-3, 5]])
    embedder = simplex2.Embedder()
    output = embedder.embed(some_input)

    print("Embedding:")
    print(f"{some_input}")
    print("leads to:")
    print(f"{output}")

    print("Decoding then gives:")

    decoder = Decoder()
    decoded_input = decoder.decode(2, 2, output[0])
    print(f"{decoded_input}")
    
        
        
    
        
        
        
            
                

        
    