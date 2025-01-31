# Specify ASCII replacements for specific Unicode characters or sequences of
# Unicode characters. This list is not intended to be exhaustive, but does
# include all the Unicode characters so far encountered in Licor data files. It
# also includes other characters that should not be allowed in column names or
# values, such as tabs.
#
# Important note about this code: the order of these replacements matters.
# For example, the replacement for `<c2><b2>` must occur after the replaement
# for `<e2><81><bb><c2><b2>`; otherwise, instances of `<e2><81><bb><c2><b2>`
# would get partially replaced by the replacement for `<c2><b2>`.
UNICODE_REPLACEMENTS <- data.frame(
    rbind(
        c('<ce><94>',             'Delta'),
        c('<ce><b1>',             'alpha'),
        c('<e2><81><bb><c2><b2>', '^(-2)'),
        c('<e2><81><bb><c2><b9>', '^(-1)'),
        c('<c2><b5>',             'micro'),
        c('<c2><b0>',             'degrees '),
        c('<c2><b2>',             '^2'),
        c('\t',                   ' ')
    ),
    stringsAsFactors = FALSE
)
colnames(UNICODE_REPLACEMENTS) <- c('pattern', 'replacement')

# Define a helping function for replacing Unicode characters.
#
# In R, some problems occur for column names or units that include Unicode
# characters such as Greek letters. We use a combination of the `iconv` and
# `gsub` functions to replace them by equivalent phrases, following a strategy
# discussed on
# https://stackoverflow.com/questions/36108790/trouble-with-strings-with-u0092-unicode-characters
replace_unicode <- function(strings)
{
    strings <- iconv(
        strings,
        '',
        'ASCII',
        'byte'
    )

    for (i in seq_along(UNICODE_REPLACEMENTS[['pattern']])) {
        strings <- gsub(
            x = strings,
            pattern = UNICODE_REPLACEMENTS[['pattern']][i],
            replacement = UNICODE_REPLACEMENTS[['replacement']][i]
        )
    }

    return(strings)
}
