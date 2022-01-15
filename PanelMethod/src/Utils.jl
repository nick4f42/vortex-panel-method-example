module Utils

export cdot, ccross, normalize

cdot(a::Number, b::Number) = real(conj(a) * b)
ccross(a::Number, b::Number) = imag(conj(a) * b)
normalize(z::Number) = z / abs(z)

end # module
