require "./random_number_generator"

module Fp

	def pow(base, exp, mod)
    ans = 1
    while exp > 0
      ans = (ans * base) % mod if (exp & 1) == 1
      exp >>= 1
      base = base**2 % mod
    end
    ans
  end

	def primitive_root(q)
		pr = 1
		f = 1
		while f == 1
			f = 0
			pr += 2
			rem = q % pr
			for i in 0..pr-1 do
				if rem == i**2 % pr
					f = 1
				end
			end
		end
		pr
	end

	module_function :pow, :primitive_root

end

class PolynomialRing

	include Fp

	def initialize(q:, n:)
		@q = q
		@n = n
		@prim_root = primitive_root(@q)
		@n_inv = pow(@n, @q-2, @q)
		@stage = Math.log2(@n).to_i
		@omega = [1, pow(@prim_root, (@q-1)/2/@n, @q)]
		for i in 2..2*@n do
			@omega += [@omega[i-1]*@omega[1] % @q]
		end
	end

	def uniform_sampling
		a = Array.new(@n)
		a.map!{ |e| RNG::Generator.rand(@q) }
	end

	def dwt(poly)
		for i in 0..@n-1 do
			poly[i] = poly[i]*@omega[i] % @q
		end
		base = @n
		for stage in 0..@stage-1 do
			for i in 0..@n/base-1 do
				offset = i*base;
				for j in 0..base/2-1 do
					radix2_a = poly[offset+j]
					radix2_b = poly[offset+j+base/2]
					poly[offset+j] = (radix2_a + radix2_b) % @q
					poly[offset+j+base/2] = (radix2_a - radix2_b) * @omega[j*2*@n/base] % @q
				end
			end
			base /= 2
		end
	end

	def idwt(poly)
		base = 2
		for stage in 0..@stage-1 do
			for i in 0..@n/base-1 do
				offset = i*base;
				for j in 0..base/2-1 do
					radix2_a = poly[offset+j]
					radix2_b = poly[offset+j+base/2] * @omega[2*@n-j*2*@n/base] % @q
					poly[offset+j] = (radix2_a + radix2_b) % @q
					poly[offset+j+base/2] = (radix2_a - radix2_b) % @q
				end
			end
			base *= 2
		end
		for i in 0..@n-1 do
			poly[i] = poly[i]*@n_inv*@omega[2*@n-i] % @q
		end
	end

	def conv(a, b)
		c = Array.new(@n, 0)
		for i in 0..@n-1 do
			for j in 0..@n-1 do
				if i + j >= @n
					c[i+j-@n] = (c[i+j-@n] - a[i]*b[j]) % @q
				else
					c[i+j] = (c[i+j] + a[i]*b[j]) % @q
				end
			end
		end
		c
	end

	def inner_prod(a, b)
		c = Array.new(@n)
		for i in 0..@n-1 do
			c[i] = a[i]*b[i] % @q
		end
		c
	end

end
