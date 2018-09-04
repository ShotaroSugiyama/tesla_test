require "./const"
require "./polynomial_ring"
require "digest/sha2"

class ChaCha

  def initialize(seed)
    @x = [0x61707876, 0x3320646e, 0x79622d32, 0x6b206574]
    @x += seed
    @x += [0, 0, 0, 0]
    @y = @x.map { |e| e }
  end

  def quarter_round(a, b, c, d)
    @x[a] = (@x[a] + @x[b]) % 2**32
    @x[d] ^= @x[a]
    @x[d] = (@x[d] << 16) | (@x[d] >> (32-16))
    @x[c] = (@x[c] + @x[d]) % 2**32
    @x[b] ^= @x[c]
    @x[b] = (@x[b] << 12) | (@x[b] >> (32-12))
    @x[a] = (@x[a] + @x[b]) % 2**32
    @x[d] ^= @x[a]
    @x[d] = (@x[d] << 8) | (@x[d] >> (32-8))
    @x[c] = (@x[c] + @x[d]) % 2**32
    @x[b] ^= @x[c]
    @x[b] = (@x[b] << 7) | (@x[b] >> (32-7))
  end

  def update
    (0..9).each do |i|
      # Column Round
      quarter_round(0, 4, 8, 12)
      quarter_round(5, 9, 13, 1)
      quarter_round(10, 14, 2, 6)
      quarter_round(15, 3, 7, 11)
      # Diagonal Round
      quarter_round(0, 5, 10, 15)
      quarter_round(1, 6, 11, 12)
      quarter_round(2, 7, 8, 13)
      quarter_round(3, 4, 9, 14)
    end
    @x.map!.with_index { |e, i| (e + @y[i]) % 2**32 }
  end

end

module TESLAFunctions

  Generator = Random.new(0)

  def uniform_sampling
		(0..Const::N-1).map { |i| Generator.rand(Const::Q) }
	end

  def bounded_sampling
    (0..Const::N-1).map { |i|  Generator.rand(-Const::B..Const::B) % Const::Q}
  end

  def gaussian_sampling
    (0..Const::N-1).map do |i|
      x = Generator.rand
      y = Generator.rand
      z = Math.sqrt(-2.0*Math.log(x))*Math.sin(2*Math::PI*y)
      (z*Const::Sigma).to_i % Const::Q
    end
  end

  def d_rounding(v)
    v.map do |e|
      tmp = e % 2**(32-Const::D)
      tmp -= 2**Const::D if tmp > 2**(Const::D-1)
      (e - tmp) >> Const::D
    end
  end

  def hash512(v1, v2, message)
    d = Digest::SHA512.hexdigest(v1.pack("c*")+v2.pack("c*")+message)
    d = d.each_char.each_slice(8).map(&:join)
    d.map { |e| e.hex }
  end

  def hash_f(sig_c)
    shift = Array.new(Const::OmegaPr)
    (0..15).each do |i|
      shift[3*i] = sig_c[i] & 0x3ff;
      shift[3*i+1] = (sig_c[i] >> 10) & 0x3ff;
      shift[3*i+2] = (sig_c[i] >> 20) & 0x3ff;
    end
    """
    (0..4).each do |i|
      shift[48] += (sig_c[i] >> (30 - 2*i)) & (0x3 << 2*i);
      shift[49] += (sig_c[i+5] >> (30 - 2*i)) & (0x3 << 2*i);
      shift[50] += (sig_c[i+10] >> (30 - 2*i)) & (0x3 << 2*i);
    end

    uint32_t count = 0;
    for (i = 0; i < params.omega_pr; i++) {
      if(fc[shift[i]] == 0) {
        fc[shift[i]] = 1;
        count++;
      }
      if(count == params.omega) {
        break;
      }
    }

    free(shift);

    if (count != params.omega) {
      return 1;
    } else {
      return 0;
    }
    """
  end

  module_function :uniform_sampling, :bounded_sampling, \
  :gaussian_sampling, :d_rounding, :hash512, :hash_f

end

class TESLA256

  class << self
    include TESLAFunctions

    def key_gen
      key = {}
      key[:s] = gaussian_sampling
      key[:e1] = gaussian_sampling
      key[:e2] = gaussian_sampling

      key[:seed] = (0..7).map { |e| Generator.rand(2**32) }
      a1 = []
      a2 = []
      chacha = ChaCha.new(key[:seed])
      (0..Const::N/16-1).each do |i|
        a1 += (chacha.update).map { |e| e % Const::Q }
        a2 += (chacha.update).map { |e| e % Const::Q }
      end

      s = key[:s].map { |e| e }
      poly_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      poly_ring.dwt!(s)
      poly_ring.dwt!(a1)
      poly_ring.dwt!(a2)
      t1 = poly_ring.inner_prod(s, a1)
      t2 = poly_ring.inner_prod(s, a2)
      poly_ring.idwt!(t1)
      poly_ring.idwt!(t2)
      key[:t1] = poly_ring.add(t1, key[:e1])
      key[:t2] = poly_ring.add(t1, key[:e2])

      key
    end

    def sign(key, message)
      key[:seed] = (0..7).map { |e| Generator.rand(2**32) }
      a1 = []
      a2 = []
      chacha = ChaCha.new(key[:seed])
      (0..Const::N/16-1).each do |i|
        a1 += (chacha.update).map { |e| e % Const::Q }
        a2 += (chacha.update).map { |e| e % Const::Q }
      end

      poly_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      poly_ring.dwt!(a1)
      poly_ring.dwt!(a2)

      while true
        r = bounded_sampling
        poly_ring.dwt!(r)

        v1 = poly_ring.inner_prod(a1, r)
        v2 = poly_ring.inner_prod(a2, r)

        break

      end
    end

    def verify(key, message)

    end
  end

end
