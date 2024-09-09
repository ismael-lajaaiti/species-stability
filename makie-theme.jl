set_theme!(theme_minimal())

# Save figure with physical dimension.
cm_to_pt = 28.3465 # 1 cm = 28.3465 pt, cf. CairoMakie documentation.
single_column_width = 8.2 # In centimeners, cf. Ecology Letter guidelines.
two_third_page_width = 11
full_page_width = 17.3
width_height_ratio = 1.5

function label_panels!(fig, nrow, ncol)
    letters = string.(collect('A':'Z')) # To label figure panels.
    for (idx, label) in zip(CartesianIndices((ncol, nrow)), letters)
        println(idx)
        i, j = idx[2], idx[1]
        layout = fig[i, j] = GridLayout()
        Label(
            layout[1, 1, TopLeft()],
            label;
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :right,
        )
    end
end

function save_figure(filename, fig, size)
    for fmt in ["pdf", "png", "eps"]
        save("$filename.$fmt", fig; size = size, pt_per_unit = 1)
    end
end
