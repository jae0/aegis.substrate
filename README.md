Basic environment for using aegis*, bio*, and stm* packages. This provides the required functions to bootstrap and interact with each package.

To install you need to bootstrap from github directly:

```
  devtools::install_github( "jae0/aegis" )
```

or you can add the following into your Rprofile:

```
pkgsInstalled = .packages(all.available = TRUE)
if (!"aegis" %in% pkgsInstalled ) {
  message( "The package, aegis is missing, please install")
  message( "Install right now? (y/n):")
  o = readline()
  if (o=="y"){
    if (!"devtools" %in% pkgsInstalled ) {
      install.packages("devtools", dependencies=TRUE, ask=FALSE)
    }
    require( devtools)
    install_github( "jae0/aegis")
  }
}
```


Then, you need to have an Rprofile set up properly such as:

```.
libPaths("~/R")
homedir = path.expand("~")
tmpdir = file.path( homedir, "tmp" )
work_root = file.path( homedir, "work" )    ### replace with correct path to work directory (local temporary storage)
code_root = file.path( homedir, "bio" )   ### replace with correct path to the parent directory of your git-projects
data_root = file.path( homedir, "bio.data" )   ### replace with correct path to your data

# store your passwords and login here and make sure they are secure
passwords = file.path( homedir, ".passwords" )
if (file.exists(passwords)) source( passwords )

require( aegis )
```


Thereafter, you can used the bootstrapped environment to install the other basic tools:

```
  aegis::project.libraryInstall()
```

If you have a local git clone of the required packages, you can install with:

```
  aegis::project.libraryInstall(local=TRUE)

```

