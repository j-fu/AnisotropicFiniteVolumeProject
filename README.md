# Finite Volumes for anisotropic problems


## Supporting information
- [Some resources on git](https://j-fu.github.io/marginalia/git/basics/)
- [Some information on working with Julia](https://j-fu.github.io/marginalia/julia)
  - This respository is organized along the lines of the [project workflow](https://j-fu.github.io/marginalia/julia/project-workflow).

## Very basic git workflow
These concern working from the command line.
On Mac or Windows, consider installing SourceTree as a GUI for git

- Cloning the git repository:
     
      $ git clone https://github.com/j-fu/anisofv

- Getting new files:
     
      $ git pull

- Adding and conveying new or modified file:
      
      $ git add <filename>
      $ git commit <filename>
      $ git push
     
## Setting up the  project
Ensure you have installed at least Julia 1.7

This installs all necessary Julia packages and  should be always run after `git pull` in the repository directory: 
 
      $ julia --project=.
      julia> using Pkg
      julia> Pkg.instantiate()


## Running a notebook in `./notebooks`

      $ julia --project=.
      julia> using Pluto
      julia> Pluto.run()

This will open a new tab in your browser.
Open `./notebooks/stationary.jl` in the file dialog.

Alternatively, you can start a notebook directly:
      
      $ julia
      julia> using Pluto
      julia> Pluto.run(notebook="notebooks/stationary.jl)

If you want to change things in a notebook, please consider copying it first to another notebook in the notebooks directory.





