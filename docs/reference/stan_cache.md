# Location of the compiled Stan model cache

Models fit with `cmd_stan = TRUE` are compiled once and the executables
are reused in later R sessions. The cache lives in the directory given
by `options(stanimal.cache_dir = )` or the environment variable
`STANIMAL_CACHE_DIR`, and defaults to the user cache directory returned
by `tools::R_user_dir("stanimal", "cache")`. Executables are specific to
the CmdStan version and platform; after upgrading CmdStan, call
`stanimal_clear_cache()`.

## Usage

``` r
stan_cache_dir()

stanimal_clear_cache()
```

## Value

`stan_cache_dir()` returns the (created) cache directory.
`stanimal_clear_cache()` deletes its contents and returns the number of
files removed, invisibly.
