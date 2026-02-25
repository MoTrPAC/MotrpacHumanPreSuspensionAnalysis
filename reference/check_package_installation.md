# Check package installation

Several secondary R packages are used and it's preferable to install
them only if they are needed. Placing this check in the body of a
function will provide helpful instructions for cases when a user
attempts to run a function but doesn't have the required suggested
package(s). (The main inspiration of this was because a lot of people
had issues installing some packages needed for only specific analysis)

## Usage

``` r
check_package_installation(pkg, fun = NULL, task = NULL)
```

## Arguments

- pkg:

  (*character*) the name of the suggested R package

- fun:

  (*character*) the name of the function that uses `pkg`

- task:

  (*character*; *optional*) the task that `fun` completes.

## Value

`TRUE`, invisibly, if `pkg` is installed.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- MotrpacHumanPreSuspensionAnalysis::check_package_installation(pkg = 'base', fun = 'mean')
result
} # }
```
