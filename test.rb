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
      p TESLA256::verify(key, signature, message)
    end

  end

end

if __FILE__ == $0
  p TESLA256::verify(TestVector::Key, TestVector::Signature, "message\x00")
end
