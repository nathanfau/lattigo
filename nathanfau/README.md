# Lattigo Fork

## Overview

This repository is a fork of Lattigo v6 made on the 4th of july 2026. 

Please note that all my modifications and my implementation are located in the `nathanfau/` directory. The rest of the original Lattigo codebase remains entirely untouched.

## Packages

The `nathanfau/` directory contains the following packages. They are designed to be completely independent of one another:

* `aes`
* `bbbts`
* `bitbatching`
* `cleaning`
* `convctx`
* `debug`
* `trigo`
* `utils`

## Testing

Packages can be tested individually. You can run the tests for any specific package using the standard Go testing tool:

```bash
go test ./nathanfau/<package_name>/ -v
```

LLMs were used occasionally to help write or refactor code, and to assist with the documentation.