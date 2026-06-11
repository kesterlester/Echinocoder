import numpy as np
from bokeh.io import curdoc
from iwpc.scalars.scalar import Scalar
from iwpc.scalars.scalar_function import ScalarFunction
from iwpc.visualise.bokeh_function_visualiser_2D import BokehFunctionVisualiser2D

def make_input_scalars(n, k):
    scalars = []
    for i in range(n):
        for j in range(k):
            scalars.append(Scalar(f"Vector {i} feature {j}", bins=np.linspace(-5, 5, 100)))
    return scalars

def bla(x, i):
    # print(x)
    return x[:, i]

def make_output_scalars(R):
    return [
        ScalarFunction(
            lambda x, i=i: bla(x, i),
            label=f"Output {i}"
        ) for i in range(R)
    ]

def evaluate_embedding(x, n, k, embedder):

    ##### Lines starting like this comment out old caching mechanism
    #####
    ##### #print("Cache test")
    ##### if evaluate_embedding.last_x is not None and len(x)==len(evaluate_embedding.last_x):
    #####     if (x==evaluate_embedding.last_x).all() and (n==evaluate_embedding.last_n) and (k==evaluate_embedding.last_k):
    #####         # We hit the cache!
    #####         #print("Cache hit")
    #####         return evaluate_embedding.last_ret
    #####     else:
    #####         #print("Cache fail")
    #####         pass

    outs = []
    for sample in x:
        sample = sample.reshape(n, k)
        embedding = np.asarray(embedder.embed(sample)[0], dtype=float)
        assert len(embedding) == embedder.size_from_n_k(n,k)
        outs.append(embedding)

    ##### evaluate_embedding.last_x = x.copy()
    ##### evaluate_embedding.last_n = n
    ##### evaluate_embedding.last_k = k
    ##### evaluate_embedding.last_ret = np.asarray(outs)
    #####
    ##### return evaluate_embedding.last_ret

    return np.asarray(outs)

def add_bokeh_root():
    n = 3
    k = 2
    
    # Uncomment one of the next few pairs of lines:
    
    # from Cinf_numpy_polynomial_encoder_for_array_of_reals_as_multiset import Encoder
    # embedder = Encoder()
    
    # from Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset import Embedder
    # embedder = Embedder()
    
    # from Historical.C0_simplicialComplex_embedder_1_for_array_of_reals_as_multiset import Embedder
    # embedder = Embedder()
    
    #from C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset import Embedder
    #embedder = Embedder()
    
    from C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset import Embedder
    embedder = Embedder()
    
    # from C0HomDeg1_conjectured_dotting_embedder_for_array_of_reals_as_multiset import Embedder
    # embedder = Embedder(n=n, k=k)
    
    evaluate_embedding.last_x=None
    evaluate_embedding.last_n=None
    evaluate_embedding.last_k=None
    evaluate_embedding.last_outs=None
    
    bokeh_vis = BokehFunctionVisualiser2D(
        lambda x: evaluate_embedding(x, n, k, embedder),
        make_input_scalars(n, k),
        make_output_scalars(embedder.size_from_n_k(n,k)),
        center_point=10 * (np.random.random(size=n * k) - 0.5),
        panel_1d_kwargs={'use_points': True},
        use_points_for_xsecs=True,
    )
    curdoc().add_root(bokeh_vis.root)

add_bokeh_root()


