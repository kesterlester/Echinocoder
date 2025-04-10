#!/usr/bin/env python
from fractions import Fraction
from itertools import pairwise
import numpy as np

def numpy_array_of_frac_to_str(arr : np.array):
    """
    Numpy arrays of fractions print in a really ugly way. This method is designed to output a human-readable
    array representation that looks more readable. Its output is only for human consumption, not archival, so
    don't try to read it back in. If the array happens not to conain fractions, it should still cope ... it will just
    print the array as usual.

    Args:
        arr: array to print

    Returns:

    """

    # This next line is heuristic but good enough for now as it will distinguish float and integer from Fraction
    # (which is all we need it to do).
    # is_array_of_fractions = ( arr.dtype == "O" )

    # if not is_array_of_fractions:
    #  Fall back to normal str
    #   return str(arr)

    ans = "["
    for row in arr:
        ans += "["
        for elt in row:
            ans += str(elt) + ", "
        ans += "],  "
    ans += "]"
    return ans

def pretty_print_lin_comb(lin_comb):
    for coeff, basis_elt in lin_comb:
        print(float(coeff), numpy_array_of_frac_to_str(basis_elt))


class EncDec:

    def __init__(self, pass_forward, pass_backward):
        self.pass_forward = pass_forward
        self.pass_backward = pass_backward

    def get_default_output_dict(self, input_dict, debug=False):
        return EncDec.get_default_dict(input_dict, self.pass_forward, debug)

    def get_default_input_dict(self, output_dict, debug=False):
        return EncDec.get_default_dict(output_dict, self.pass_backward, debug)


    @staticmethod
    def get_default_dict(some_dict, some_pass_through, debug=False):

        if some_pass_through is False or some_pass_through is None or not some_pass_through:
            ret_dict = dict()
            if debug:
                print(f"some_pass_through was False or None")
        elif some_pass_through is True:
            ret_dict = some_dict
            if debug:
                print(f"some_pass_through was True")
        else:
            # Assime some_pass_through is a list
            if debug:
                print(f"some_pass_through is assumed to be a list")
            ret_dict = {k: v for k, v in some_dict.items() if k in some_pass_through}

        if debug:
            print(f"ret_dict initial value is\n{ret_dict}")

        return ret_dict


    def encode(self, input_dict):
        """
        Encoder-Decoders should implement.
        input_dict is a dictionary of key,value pairs.
        output_dict should be a dictionary of key,value pairs.
        encode maps input -> output
        """
        raise NotImplementedError()

    def decode(self, output_dict):
        """
        Encoder-Decoders should implement.
        input_dict is a dictionary of key,value pairs.
        output_dict should be a dictionary of key,value pairs.
        decode maps output -> input
        """
        raise NotImplementedError()

class PassThrough(EncDec):

    def __init__(self):
        super().__init__(pass_forward=False, pass_backward=False)

    def encode(self, input_dict):
        return input_dict

    def decode(self, output_dict):
        return output_dict

class MergeLinCombs(EncDec):
    def __init__(self, input_lin_comb_names, output_lin_comb_name, pass_forward=None, pass_backward=None):
        """
        Merges linear combinations in the input dictionary and writes them to the output dictionary.
        E.g. could combing  input["diff1"]=[(2,"i"), (5,"j")] and input["offset"] = [(3,"k")] into
        output["out"] = [(2,"i"), (5,"j"), (3,"k")].

        Args:
            input_lin_comb_names: a list of lin_comb names, e.g. ["diff1", "offset"].
            output_lin_comb_name: the name of the output lin comb, e.g. "out".
        """
        super().__init__(pass_forward, pass_backward)
        self.input_lin_comb_names = input_lin_comb_names
        self.output_lin_comb_name = output_lin_comb_name

    def encode(self, input_dict, debug=False):
        out_dict = self.get_default_output_dict(input_dict, debug)
        out_lin_comb = []
        for lin_comb_name in self.input_lin_comb_names:
            lin_comb = input_dict[lin_comb_name]
            out_lin_comb.extend(lin_comb)
        out_dict[self.output_lin_comb_name] = out_lin_comb
        return out_dict


class ArrayToLinComb(EncDec):
    def __init__(self, input_array_name, output_lin_comb_name, pass_forward=None, pass_backward=None):
        super().__init__(pass_forward, pass_backward)
        self.input_array_name = input_array_name
        self.output_lin_comb_name = output_lin_comb_name

    def encode(self, input_dict, debug=False):
        if self.input_array_name not in input_dict:
            raise ValueError(f"Expected {self.input_array_name} in {input_dict}.")

        output_dict = self.get_default_output_dict(input_dict, debug)
        arr = np.asarray(input_dict[self.input_array_name])
        result = []
        for index, value in np.ndenumerate(arr):
            basis = np.zeros_like(arr)
            basis[index] = 1
            result.append((value, basis))

        output_dict[self.output_lin_comb_name] = result
        if debug:
            print(f"About to return {output_dict}")
        return output_dict

class Chain(EncDec):
    def __init__(self, encoder_decoder_list):
        self.encoder_decoder_list = encoder_decoder_list

    def encode(self, input_dict, debug=False):
        for encoder in self.encoder_decoder_list:
            if debug:
                print(f"Chain is about to run {encoder.__class__}")
            input_dict = encoder.encode(input_dict, debug=debug)
        return input_dict

    def decode(self, output_dict, debug=False):
        for encoder in reversed(self.encoder_decoder_list):
            output_dict = encoder.encode(output_dict, debug=debug)
        return output_dict


class BarycentricSubdivide(EncDec):
    """
    Encoding:
        * looks for a dictionary entry named self.input_name which is assumed to be a point in barycentric coordinates,
          i.e. a linear combination of basis element representing the vertices of a simplex,
        * calculates how this poing would be expressed ith respect to a different basis corresponding to a barycentric
          subdivision of the simplex,
        * creates a dictionary with an entry for self.output_name carrying the resulting linear combination, and
        * if "pass_through" is true, the encoder (and decode) will pass through all other elements of the input
          dictionary which do not result in overwriting of one of lin comb output.
        * If preserve_scale is True (default) then the sum of the coeffiencients is preserved. Equivalently, the one
          norm of each basis vector iw preserved at 1 if already at 1.

        A linear combination is assumed to be represented as a list of pairs -- with the first element of each pair
        being the coefficient and the second element of each pair being the corresponding basis verctor. I.e.
        [ (2, ei), (1, ej), (3, ek)]
        might represent the vector (2,1,3) with respect to the usual cartesian basis.

    Decoding:
        Reverse of encoding.
    """
    def __init__(self, input_name, diff_output_name, offset_output_name, pass_forward=None, pass_backward=None, preserve_scale=True):
        """

        Args:
            input_name: name of lin-comb to encode to on
            diff_output_name: name of lin-comb into which to encode everything EXCEPT the constant offset term.
                               This lin-comb has only non-negative coefficients.
            offset_output_name: name of the lin-comb to encode the constant offset term.
                                This lin-comb may have coefficients of any sign.
            pass_through: whether or not to pass metadata in the input to the output (or vice versa).
                          If True, then pass everything through.
                          If False, then pass nothing through.
                          If list, then pass through dict items whose names are in the list
        """
        
        super().__init__(pass_forward, pass_backward)
        self.input_name = input_name
        self.diff_output_name = diff_output_name
        self.offset_output_name = offset_output_name
        self.common_output_name = diff_output_name if diff_output_name==offset_output_name else None
        self.preserve_scale = preserve_scale


    def encode(self, input_dict, pass_through=None, debug=False):
        if debug:
            print(f"input_dict is\n{input_dict}")

        output_dict = self.get_default_output_dict(input_dict, debug)

        if not self.input_name in input_dict:
            raise ValueError(f"Expected {self.input_name} in {input_dict}.")

        input_lin_comb = input_dict[self.input_name]
        """
        A linear combination is assumed to be represented as a list of pairs -- with the first element of each pair
        being the coefficient and the second element of each pair being the corresponding basis vector. I.e.
        [ (2, "i"), (1, "j"), (3, "k")]
        might represent the vector (2,1,3) with respect to the usual cartesian basis.
        """
        if debug:
            print(f"input_lin_comb is\n{input_lin_comb}")

        if not input_lin_comb:
            # List for p is empty, so
            if self.common_output_name:
                output_dict[self.common_output_name] = list()
            else:
                output_dict[self.diff_output_name] = list()
                output_dict[self.offset_output_name] = list()

            if debug:
                print(f"About to return \n{output_dict}")
            return output_dict

        assert len(input_lin_comb) >= 1

        # Sort by coefficient in linear combination, big -> small
        # We need a list below as we will consume it twice when generating the diff_lin_comb
        sorted_lin_comb = sorted(input_lin_comb, key=lambda x: x[0], reverse=True)
        if debug:
            print(f"sorted_lin_comb is\n{sorted_lin_comb}")

        coeffs = [ x for x, _ in sorted_lin_comb ]
        basis_vecs = [x for _ , x in sorted_lin_comb]

        """
        class ZeroObject(object):
            def __add__(self, other):
                return other
        """

        if self.preserve_scale:
            diff_lin_comb = list(((fac := (i+1))*(x-y), sum(basis_vecs[:i+1], start=0*basis_vecs[0]+Fraction())/fac) for i, (x,y) in enumerate(pairwise(coeffs)))
            offset_lin_comb = [((fac := len(basis_vecs))*coeffs[-1], sum(basis_vecs, start=0*basis_vecs[0]+Fraction())/fac)]
        else:
            diff_lin_comb = list(((x-y), sum(basis_vecs[:i+1], start=0*basis_vecs[0])) for i, (x,y) in enumerate(pairwise(coeffs)))
            offset_lin_comb = [(coeffs[-1], sum(basis_vecs, start=0*basis_vecs[0]))]

        if debug:
            print(f"diff_lin_comb is\n{diff_lin_comb}")
            print(f"offset_lin_comb is\n{offset_lin_comb}")

        if self.common_output_name:
            output_dict[self.common_output_name] = diff_lin_comb + offset_lin_comb
        else:
            output_dict[self.diff_output_name] = diff_lin_comb
            output_dict[self.offset_output_name] = offset_lin_comb
        if debug:
            print(f"About to return \n{output_dict}")

        return output_dict

def tost():
    print("###########################################")

    pass_through = PassThrough()

    input = {"hello": "world"}
    enc = pass_through.encode(input)
    dec = pass_through.decode(enc)
    assert dec == input
    assert enc == input

    print("###########################################")

    subdivide = BarycentricSubdivide("p","p", "q")

    input_dict = {"p": [(-1, np.array([1, 0, 0])),
                        (-7, np.array([0, 1, 0])),
                        (10, np.array([0, 0, 1]))],
                  "metadata1": "moo1",
                  "metadata2": "moo2",
                  "metadata3": "moo3",
                  }
    enc = subdivide.encode(input_dict, debug=True)


    print("###########################################")

    input_dict = { "arr" : np.asarray([[ 4, 2],
                                       [-3, 5],
                                       [ 8, 9],
                                       [ 2 ,7]]) }
    array_to_lin_comb = ArrayToLinComb("arr", "lin_comb" )

    enc = array_to_lin_comb.encode(input_dict)
    print(f"=======================\nArray to lin comb made encoded")
    print(input_dict)
    print("to")
    print(f"{enc}")

    print("###########################################")
    simplex1_bit = Chain([
        BarycentricSubdivide("set", "first_diffs", "offset"),
        BarycentricSubdivide("first_diffs", "second_diffs", "second_diffs", pass_forward="offset")
    ])

    input_dict = {
                  "set": [(-1, np.array([1, 0, 0])),
                          (-7, np.array([0, 1, 0])),
                          (10, np.array([0, 0, 1]))],
                  "metadata1": "moo1",
                  "metadata2": "moo2",
                  "metadata3": "moo3",
                  }
    enc = simplex1_bit.encode(input_dict, debug=True)
    print(f"=======================\nSimplex1 as a chain encoded")
    print(input_dict)
    print("to")
    print(f"{enc}")



    print("###########################################")
    simplex1_different_bit = Chain([
        ArrayToLinComb(input_array_name="set", output_lin_comb_name="lin_comb_0"),
        BarycentricSubdivide("lin_comb_0", "lin_comb_1_first_diffs", "offset", preserve_scale=False),
        BarycentricSubdivide("lin_comb_1_first_diffs", "lin_comb_2_second_diffs",
                             "lin_comb_2_second_diffs", pass_forward="offset", pass_backward="offset", preserve_scale=False),
        MergeLinCombs(["lin_comb_2_second_diffs", "offset"], "lin_comb_3"),
    ])

    input_dict = {
                  "set" : np.asarray([[ 4, 2],
                                      [-3, 5],
                                      [ 8, 9],
                                      [ 2 ,7]]),
                  "metadata1": "moo1",
                  "metadata2": "moo2",
                  "metadata3": "moo3",
                  }
    enc = simplex1_different_bit.encode(input_dict, debug=True)
    print(f"=======================\nSimplex1 as a chain encoded")
    print(numpy_array_of_frac_to_str(input_dict["set"]))
    print("to")
    #print(f"{enc}")

    lin_comb_3 = enc["lin_comb_3"]
    #print(f"Note that lin_comb_3 is")
    pretty_print_lin_comb(lin_comb_3)
    print("and the (non-offset) differences are")
    [ print(numpy_array_of_frac_to_str(tmp:=b-a), " with ", np.sum(tmp)," ones in it") for a,b in list(pairwise( [a for _,a in lin_comb_3 ]))[:-1] ]


    print("###########################################")

    print("All tests passed.")


if __name__ == "__main__":
    tost()
