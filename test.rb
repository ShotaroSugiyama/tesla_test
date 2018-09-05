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
      p_ring = PolynomialRing.new(q: Const::Q, n: Const::N)

      e1 = key[:e1]
      print10(e1)
      e2 = p_ring.dwt(e1)
      print10(e2)
      e1 = p_ring.idwt(e2)
      print10(e2)
    end

  end

end

if __FILE__ == $0
  key = TESLA256::key_gen
  message = "message\x00"
  signature = TESLA256::sign(key, message)
  p TESLA256::verify(key, signature, message)
end
