"""
    phase_diff(field, dir)

Compute point-to-point difference of phase on a grid along a direction

Returns an array of size `(size(field)[1], size(field)[2], size(field)[2])`
containing gradients in direction `dir`.
"""
function phase_diff(field, dir)
    out = similar(field, Real)

    out .= angle.(circshift(field, dir) ./  field)

    out
end

"""
    phase_diff(field::PencilArray, dir)

Compute point-to-point difference of phase on a grid along a direction

Returns an array of size `(size(field)[1], size(field)[2], size(field)[2])`
containing gradients in direction `dir`.
"""
function phase_diff(field::PencilArray, dir)
    out = similar(field, Real)

    # TODO diff at boundaries
    out .= 0
    out[1:end-dir[1], 1:end-dir[2], 1:end-dir[3]] .= angle.(circshift(field[1:end-dir[1], 1:end-dir[2], 1:end-dir[3]], dir) ./ field[1+dir[1]:end, 1+dir[2]:end, 1+dir[3]:end])

    out
end

"""
    normed_max_phase_grad(grids)

Compute maximum phase gradient of a grid

Normalised to ignore large gradients in regions with low density.  These tend
to be anomalous.
"""
function max_normed_phase_grad(grids)
    DENSITY_THRESHOLD = 1e-6
    tmp = similar(grids.ψx, Real)
    n = ndims(tmp)

    max_grads = zeros(n)

    shift = zeros(Int, n)
    shift[1] = 1

    for i in 1:n
        dir = circshift(shift, i)
        tmp .= abs.(phase_diff(grids.ψx, dir))
        tmp[grids.ρx / maximum(grids.ρx) .< DENSITY_THRESHOLD] .= NaN
        max_grads[i] = maximum(tmp)
    end

    maximum(max_grads)

end
