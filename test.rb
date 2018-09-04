require "./const"
require "./polynomial_ring"
require "./tesla"

class Test

	class << self

		def print10(poly)
			for i in 0..9 do
				puts "#{i} : #{poly[i]}"
			end
			puts ""
		end

		def test_dwt
			key = TESLA256.key_gen
			poly_ring = PolynomialRing.new(q: Const::Q, n: Const::N)

			e1 = key[:e1]
			e2 = key[:e2]
			e3 = poly_ring.conv(e1, e2)
			Test::print10(e3)
			poly_ring.dwt!(e1)
			poly_ring.dwt!(e2)
			e3 = poly_ring.inner_prod(e1, e2)
			poly_ring.idwt!(e3)
			Test::print10(e3)
		end

	end

end

if __FILE__ == $0
	Test.test_dwt
end
