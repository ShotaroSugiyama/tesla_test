require "./const"
require "./polynomial_ring"
require "./random_number_generator"

class Test

	class << self

		def print10(poly)
			for i in 0..9 do
				puts "#{i} : #{poly[i]}"
			end
			puts ""
		end

	end

end

if __FILE__ == $0

	Poly = PolynomialRing.new(q: Const::Q, n: Const::N)
	a = Poly.uniform_sampling
	b = Poly.uniform_sampling
	c = Poly.conv(a, b)
	Test::print10(c)
	Poly.dwt(a)
	Poly.dwt(b)
	d = Poly.inner_prod(a, b)
	Poly.idwt(d)
	Test::print10(d)
end
