# Custom data directory
Custom stellar data and model database files can always be used by providing a full path in the argument. However, users may optionally specify their own data directory using the environment variable `STARFIT_DATA` for convenience:
```shell
export STARFIT_DATA='/your/custom/data'
```
Files found in the custom data directory will take precedence over the default data directory.
Your custom data directory must have the same structure as `src/starfit/data`, i.e. it should contain the `db`, `ref`, and `stars` directories:
```shell
❯ ls
db
ref
stars
```
