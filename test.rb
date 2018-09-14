require "./const"
require "./polynomial_ring"
require "./tesla"
require "./out"

class Test

  class << self

    def sign_verify
      key = TESLA256::key_gen
      message = "message\x00"
      signature = TESLA256::sign(key, message)
      # 署名検証が通れば1、失敗した場合0
      p TESLA256::verify(key, signature, message)
    end

    def fm
      key = TestVector::Key
      signature = TestVector::Signature
      p TESLA256::verify(key, signature, "message\x00")
    end

  end

end

if __FILE__ == $0
  Test::sign_verify
  Test::fm
end
