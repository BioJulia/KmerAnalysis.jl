# Package installation

MerCounting is a component of the GenomeGraphs framework for graph based genome
assembly and analysis. So if you've installed GenomeGraphs, you should already
have it.

However, MerCounting has been designed to be generally useful in other settings
and projects too. So it can be installed on its own as well.

MerCounting is made available to install through BioJulia's package registry.

Julia by default only watches the "General" package registry, so before you
start, you should add the BioJulia package registry.

Start a julia terminal, hit the `]` key to enter pkg mode (you should see the
prompt change from `julia>` to `pkg>`), then enter the following command:

```
pkg> registry add https://github.com/BioJulia/BioJuliaRegistry.git
```

After you've added the registry, you can install GenomeGraphs from the julia REPL.
Press `]` to enter pkg mode again, and enter the following:

```
pkg> add MerCounting
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.