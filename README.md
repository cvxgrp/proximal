# Proximal operators

This "library" contains sample implementations of various proximal operators in
Matlab. These implementations are intended to be pedagogical, not the most
performant.

This code is associated with the paper 
*[Proximal Algorithms](http://www.stanford.edu/~boyd/papers/prox_algs.html)* 
by Neal Parikh and Stephen Boyd.

## Requirements

The C functions rely on the GNU Scientific Library (GSL). Some of these
functions also contain OpenMP directives to parallelize some `for` loops, so
compiling with OpenMP is optional, but some of the functions will be
substantially faster if it is used.

The Matlab function `prox_cvx.m` requires [CVX](http://cvxr.com/cvx).

## Examples

Evaluating the proximal operator of the *l1* norm via CVX and the function here:

```matlab
>> n = 100;
>> lambda = 1;
>> 
>> v = randn(n,1);
>> 
>> % CVX baseline
>> cvx_begin quiet
>>   variable x(n)
>>   minimize(norm(x,1) + (1/(2*lambda))*sum_square(x - v))
>> cvx_end
>> 
>> % Custom method
>> x2 = prox_l1(v, lambda);
>> 
>> % Comparison
>> norm(x - x2)
ans =
7.7871e-05
```

Evaluating the proximal operator of the nuclear norm:

```matlab
>> m = 10;
>> n = 30;
>> lambda = 1;
>> 
>> V = randn(m,n);
>> 
>> % CVX baseline
>> cvx_begin quiet
>>   variable X(m,n)
>>   minimize(norm_nuc(X) + (1/(2*lambda))*square_pos(norm(X - V,'fro')))
>> cvx_end
>> 
>> % Custom method
>> X2 = prox_matrix(V, lambda, @prox_l1);
>> 
>> % Comparison
>> norm(X - X2)
ans =
1.9174e-05
```

This second example shows a case where one of the arguments is a function
handle to another proximal operator.

The other Matlab functions work similarly; just use `help` in Matlab.

For a C example, see the file `example.c` in the C source directory.

## Proximal operators

The Matlab functions include the following examples:

* Projection onto an affine set
* Projection onto a box
* Projection onto the consensus set (averaging)
* Projection onto the exponential cone
* Projection onto the nonnegative orthant
* Projection onto the second-order cone
* Projection onto the semidefinite cone
* Proximal operator of a generic function (via CVX)
* Proximal operator of the *l1* norm
* Proximal operator of the max function
* Proximal operator of a quadratic function
* Proximal operator of a generic scalar function (vectorized)
* Proximal operator of an orthogonally invariant matrix function
* Precomposition of a proximal operator

## Authors

* [Neal Parikh](http://cs.stanford.edu/~npparikh)
* [Eric Chu](http://www.stanford.edu/~echu508)
* [Stephen Boyd](http://www.stanford.edu/~boyd)

## Other libraries

There are other libraries with implementations of proximal or projection
operators that may be preferable or contain more examples:

* [TFOCS](http://cvxr.com/tfocs/functions/) (see prox/proj sections) 
by S. Becker, E. Candès, and M. Grant
* [Python proximal operators](https://github.com/svaiter/pyprox) by S. Vaiter

## License

This code is released under a BSD license; see the "LICENSE" file.
