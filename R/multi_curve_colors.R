multi_curve_colors <- function() {
    c(
        "#000000",
        "#66C2A5",
        "#FC8D62",
        "#8DA0CB",
        "#E78AC3",
        "#A6D854",
        "#FFD92F",
        "#E5C494",
        "#B3B3B3",
        "#A6CEE3",
        "#1F78B4",
        "#B2DF8A",
        "#33A02C",
        "#FB9A99",
        "#E31A1C",
        "#FDBF6F",
        "#FF7F00",
        "#CAB2D6",
        "#6A3D9A",
        "#B15928",
        "#1B9E77",
        "#D95F02",
        "#7570B3",
        "#E7298A",
        "#66A61E",
        "#E6AB02",
        "#A6761D",
        "#666666"
    )
}

multi_curve_point_colors <- function() {
    c(
        multi_curve_colors()[1],
        paste0(multi_curve_colors()[seq(2, length(multi_curve_colors()))], '00')
    )
}

multi_curve_line_colors <- function() {
    lc <- multi_curve_colors()
    lc[1] <- paste0(lc[1], '00')
    lc
}

# Set point colors so that only the first curve has points; the others will all
# have points with alpha = 0 (which is fully transparent)
point_colors <- c(
    multi_curve_colors()[1],
    paste0(multi_curve_colors()[seq(2, length(multi_curve_colors()))], '00')
)

# Set line colors so that the first curve has a line, but the others don't; the
# first line will have alpha = 0 (which is fully transparent)
line_colors <- multi_curve_colors()
line_colors[1] <- paste0(line_colors[1], '00')
