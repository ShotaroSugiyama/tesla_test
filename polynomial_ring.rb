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
			(0..pr-1).map { |i| f = 1 if i**2 == q % pr}
		end
		pr
	end

	def mod_mult_tesla256(a, b)
		q = 2**28 - 2**16 + 1
		a0 = a & 0x3fff
		a1 = a >> 14
		b0 = b & 0x3fff;
		b1 = b >> 14;
		a0b0 = a0*b0;
		a1b1 = a1*b1;
		middle = (a0 + a1)*(b0 + b1) % q
		middle = (middle - a0b0 - a1b1) % q
		a0 = (a0b0 + ((middle & 0x3fff) << 14) + (middle >> 14)*0xffff) % q
		a0 = (a0 - a1b1) % q
		a1 = (((a1b1 >> 12)*0xffff % q) + ((a1b1 & 0xfff) << 16) % q) % q
		(a0 + a1) % q
	end

	module_function :pow, :primitive_root, :mod_mult_tesla256

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
		(2..2*@n).each do |i|
			@omega += [@omega[i-1]*@omega[1] % @q]
		end
	end

	def dwt!(poly)
		poly.map!.with_index { |e, i| e*@omega[i] % @q}
		base = @n
		(0..@stage-1).each do |stage|
			(0..@n/base-1).each do |i|
				offset = i*base;
				(0..base/2-1).each do |j|
					radix2_a = poly[offset+j]
					radix2_b = poly[offset+j+base/2]
					poly[offset+j] = (radix2_a + radix2_b) % @q
					poly[offset+j+base/2] = (radix2_a - radix2_b) * @omega[j*2*@n/base] % @q
				end
			end
			base /= 2
		end
	end

	def idwt!(poly)
		base = 2
		(0..@stage-1).each do |stage|
			(0..@n/base-1).each do |i|
				offset = i*base;
				(0..base/2-1).each do |j|
					radix2_a = poly[offset+j]
					radix2_b = poly[offset+j+base/2] * @omega[2*@n-j*2*@n/base] % @q
					poly[offset+j] = (radix2_a + radix2_b) % @q
					poly[offset+j+base/2] = (radix2_a - radix2_b) % @q
				end
			end
			base *= 2
		end
		poly.map!.with_index { |e, i| e*@n_inv*@omega[2*@n-i] % @q}
	end

	def conv(a, b)
		c = Array.new(@n, 0)
		(0..@n-1).to_a.repeated_permutation(2) do |i, j|
			if i + j >= @n
				c[i+j-@n] = (c[i+j-@n] - a[i]*b[j]) % @q
			else
				c[i+j] = (c[i+j] + a[i]*b[j]) % @q
			end
		end
		c
	end

  def add(a, b)
		a.map.with_index { |e, i| (e + b[i]) % @q}
  end

  def sub(a, b)
		a.map.with_index { |e, i| (e - b[i]) % @q}
  end

	def inner_prod(a, b)
		a.map.with_index { |e, i| (e * b[i]) % @q}
	end

end
