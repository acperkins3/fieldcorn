# fieldcorn
Functions useful to the fieldcorn lab

## Installation

```
devtools:install_github("acperkins3/fieldcorn")
```

## `make.IBS.table`

A function to calculate identify by state (IBS) between two or more lines from genotypes in [Hapmap format](https://bitbucket.org/tasseladmin/tassel-5-source/wiki/UserManual/Load/Load).

The current version has no issues I'm aware of but is slow. I will try vectorizing it at some point for speed.

