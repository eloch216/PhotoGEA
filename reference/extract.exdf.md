# Access or modify exdf elements

Returns or sets the values of elements in an `exdf` object.

## Usage

``` r
# S3 method for class 'exdf'
x[i, j, return_exdf = FALSE]

  # S3 method for class 'exdf'
x[i, j] <- value
```

## Arguments

- x:

  An `exdf` object.

- i, j:

  Indices specifying elements to extract or replace. Indices are
  `numeric` or `character` vectors or empty (missing) or `NULL`.

- return_exdf:

  A logical value indicating whether the return value should be an
  `exdf` object.

- value:

  Typically an array-like R object of a similar class as x.

## Details

Since an `exdf` object is actually a list of named elements, those
elements can be accessed using the `[[` or `$` operators, and a list of
all named elements can be obtained by calling `names`.

Elements of the `main_data` data frame of an `exdf` object can be
accessed and set using the `[` and `[<-` operators. When applied to an
`exdf` object, these operators are essentially shortcuts to calling the
same operators on the object's `main_data` data frame.

To create a new `exdf` object with a subset of the data contained in
another `exdf` object, the `[` operator with `return_exdf = TRUE` can be
used.

## Value

When `return_exdf` is `FALSE`, the access operator will return either a
vector or a data frame, depending on the dimension of `j`. When
`return_exdf` is `TRUE`, the access operator will return an
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) object.

## See also

[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md)

## Examples

``` r
# Create a small exdf object that includes an extra element in addition to the
# required ones (`main_data`, `units`, and `categories`).
small_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2'),
  extra_exdf_element = "This is an example of an extra exdf element"
)

# Accessing elements of `small_exdf`
names(small_exdf)     # Get the names of all elements of small_exdf
#> [1] "main_data"          "units"              "categories"        
#> [4] "extra_exdf_element"
small_exdf[['units']] # View the units using the `[[` operator
#>   A B
#> 1 m s
small_exdf$categories # View the categories using the `$` operator
#>      A    B
#> 1 Cat1 Cat2

# Accessing elements of `small_exdf$main_data`
small_exdf[,1]   # Access the first column
#> [1] 3 2 7 9
small_exdf[1,]   # Access the first row
#>   A B
#> 1 3 4
small_exdf[,'B'] # Access the column named 'B'
#> [1] 4 5 1 8
small_exdf[1,2]  # Access element 1 of column 2
#> [1] 4

# Equivalent (but longer) commands for accessing elements of `small_exdf$main_data`
small_exdf$main_data[,1]   # Access the first column
#> [1] 3 2 7 9
small_exdf$main_data[1,]   # Access the first row
#>   A B
#> 1 3 4
small_exdf$main_data[,'B'] # Access the column named 'B'
#> [1] 4 5 1 8
small_exdf$main_data[1,2]  # Access element 1 of column 2
#> [1] 4

# Replacing elements of `small_exdf$main_data`
small_exdf[,'A'] <- seq_len(4)             # Replace column A with new values
small_exdf[small_exdf[,'A'] > 2, 'B'] <- 0 # Replace some rows of column B with new values

# Creating a new exdf object with a subset of the data from small_exdf. Here we
# specify `return_exdf = TRUE` so that the `[` operator returns an exdf object
# instead of a data frame
new_exdf <- small_exdf[small_exdf[,'A'] > 2, , TRUE]
names(new_exdf) # Check that the `extra_exdf_element` is still present
#> [1] "main_data"          "units"              "categories"        
#> [4] "extra_exdf_element"
print(new_exdf) # Check that only the rows with A > 2 are included
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>   A [Cat1] (m) B [Cat2] (s)
#> 3            3            0
#> 4            4            0
```
