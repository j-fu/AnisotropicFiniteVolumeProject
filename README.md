How to use the code
===================

There are two main codes now:

notebooks/stationary.jl
notebooks/transient.jl


These are pluto notebooks which can be run in the browser in the following way:
You need to have installed Julia 1.7



Go to the folder where you checked out the code

```
$ julia --project=.
AnisoFV> using Pkg
AnisoFV> Pkg.instantiate() # this will precompile things
AnisoFV> using Pluto
AnisoFV> Pluto.run(notebook="notebooks/stationary.jl")
```


