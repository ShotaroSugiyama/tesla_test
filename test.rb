require "./const"
require "./polynomial_ring"
require "./tesla"
require "./out"

class Test

  class << self

    def print10(poly)
      for i in 0..9 do
        puts "#{i} : #{poly[i]}"
      end
      puts ""
    end

    def dwt
      p_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      x = TESLAFunctions::uniform_sampling
      print10(x)
      print10(p_ring.idwt(p_ring.dwt(x)))
    end

    def sign_verify
      key = TESLA256::key_gen
      message = "message\x00"
      signature = TESLA256::sign(key, message)
      p TESLA256::verify(key, signature, message)
    end

    def chacha
      TESLA256::key_gen
    end

    def key_check
      p_ring = PolynomialRing.new(q: Const::Q, n: Const::N)
      a1 = Array.new(Const::N, 0)
      a2 = Array.new(Const::N, 0)
      s = p_ring.dwt(TestVector::Key[:s])
      e1 = p_ring.dwt(TestVector::Key[:e1])
      e2 = p_ring.dwt(TestVector::Key[:e2])
      t1 = p_ring.dwt(TestVector::Key[:t1])
      t2 = p_ring.dwt(TestVector::Key[:t2])
      p t1 == p_ring.add(p_ring.inner_prod(a1, s), e1)
      p t2 == p_ring.add(p_ring.inner_prod(a2, s), e2)
    end

  end

end

if __FILE__ == $0
  p TESLA256::verify(TestVector::Key, TestVector::Signature, "message\x00")
end
