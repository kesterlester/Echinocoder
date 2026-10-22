"""
The ndenumerate_slice functions are supposed to work like numpu.ndenumerate, but are able to restrict themselves to part of an array 
that identified by a slice specification.
You could achieve the same effect by using ndenumerate and then filtering the output to include only things in the relevant slice,
however that would be wasteful as the rest of the array would be traversed only to be discarded. This implementation avoids that waste.
"""

import numpy as np

def ndenumerate_fixed_axis_slice(arr, axis, fixed_index):
    """
    Iterate over a slice of `arr` where one axis is fixed at `fixed_index`.
    Yields (full_index_tuple, value) pairs.
    
    Parameters:
        arr (np.ndarray): Input N-dimensional array.
        axis (int): Axis to fix.
        fixed_index (int): Index to fix on the specified axis.
        
    Yields:
        Tuple[index_tuple, value]
    """
    return ndenumerate_fixed_axes_slice(arr, {axis: fixed_index})


def ndenumerate_fixed_axes_slice(arr, fixed_axes):
    """
    Iterate over a slice of `arr` where one or more axes are fixed at given indices.
    Yields (full_index_tuple, value) pairs.
    
    Parameters:
        arr (np.ndarray): Input N-dimensional array.
        fixed_axes (dict): Dictionary mapping axis index -> fixed index.
        
    Yields:
        Tuple[index_tuple, value]
    """
    # Build the slicer with fixed axes
    slicer = []
    for axis in range(arr.ndim):
        if axis in fixed_axes:
            slicer.append(fixed_axes[axis])
        else:
            slicer.append(slice(None))
    
    # Get the sliced view
    sliced_view = arr[tuple(slicer)]

    # Iterate over the view and reconstruct full indices
    for idx, val in np.ndenumerate(sliced_view):
        full_idx = []
        dim_idx = 0
        for axis in range(arr.ndim):
            if axis in fixed_axes:
                full_idx.append(fixed_axes[axis])
            else:
                full_idx.append(idx[dim_idx])
                dim_idx += 1
        yield tuple(full_idx), val


# --- Demo ---
# --- Demo should produce output similar to this:
def demo():
    """
    Normal ndenumerate:
    (0, 0, 0) 0
    (0, 0, 1) 1
    (0, 0, 2) 2
    (0, 0, 3) 3
    (0, 1, 0) 4
    (0, 1, 1) 5
    (0, 1, 2) 6
    (0, 1, 3) 7
    (0, 2, 0) 8
    (0, 2, 1) 9
    (0, 2, 2) 10
    (0, 2, 3) 11
    (1, 0, 0) 12
    (1, 0, 1) 13
    (1, 0, 2) 14
    (1, 0, 3) 15
    (1, 1, 0) 16
    (1, 1, 1) 17
    (1, 1, 2) 18
    (1, 1, 3) 19
    (1, 2, 0) 20
    (1, 2, 1) 21
    (1, 2, 2) 22
    (1, 2, 3) 23
    
    Iterating over arr[1, :, :]:
    (1, 0, 0) 12
    (1, 0, 1) 13
    (1, 0, 2) 14
    (1, 0, 3) 15
    (1, 1, 0) 16
    (1, 1, 1) 17
    (1, 1, 2) 18
    (1, 1, 3) 19
    (1, 2, 0) 20
    (1, 2, 1) 21
    (1, 2, 2) 22
    (1, 2, 3) 23
    
    Iterating over arr[1, :, 2]:
    (1, 0, 2) 14
    (1, 1, 2) 18
    (1, 2, 2) 22
    """
    arr = np.arange(2*3*4).reshape(2, 3, 4)
    print("Normal ndenumerate:")
    for idx, val in np.ndenumerate(arr):
        print(idx, val)

    print("\nIterating over arr[1, :, :]:")
    for idx, val in ndenumerate_fixed_axis_slice(arr, axis=0, fixed_index=1):
        print(idx, val)

    print("\nIterating over arr[1, :, 2]:")
    for idx, val in ndenumerate_fixed_axes_slice(arr, fixed_axes={0: 1, 2: 2}):
        print(idx, val)


if __name__ == "__main__":
    demo()

