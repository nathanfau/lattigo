// AES-128 oracle
//
// FIPS-197 test vector:
//
//	key        = 2b7e151628aed2a6abf7158809cf4f3c
//	plaintext  = 3243f6a8885a308d313198a2e0370734
//	ciphertext = 3925841d02dc09fbdc118597196a0b32
package aes

func xtime(b byte) byte {
	if b&0x80 != 0 {
		return (b << 1) ^ 0x1b
	}
	return b << 1
}

// gfMul multiplies a and b in GF(2^8) (AES field).
func gfMul(a, b byte) byte {
	var p byte
	for i := 0; i < 8; i++ {
		if b&1 != 0 {
			p ^= a
		}
		hi := a & 0x80
		a <<= 1
		if hi != 0 {
			a ^= 0x1b
		}
		b >>= 1
	}
	return p
}

// gfInv returns the multiplicative inverse of a in GF(2^8) (0 maps to 0).
func gfInv(a byte) byte {
	if a == 0 {
		return 0
	}
	r := byte(1)
	for e := 0; e < 254; e++ {
		r = gfMul(r, a)
	}
	return r
}

var sbox [256]byte
var invSbox [256]byte

func rotl8(x byte, n int) byte { return (x << n) | (x >> (8 - n)) }

func init() {
	for i := 0; i < 256; i++ {
		q := gfInv(byte(i))
		s := q ^ rotl8(q, 1) ^ rotl8(q, 2) ^ rotl8(q, 3) ^ rotl8(q, 4) ^ 0x63
		sbox[i] = s
		invSbox[s] = byte(i)
	}
}

// Encryption round functions.

func AddRoundKey(s, rk []byte) {
	for i := range s {
		s[i] ^= rk[i]
	}
}

func SubBytes(s []byte) {
	for i := range s {
		s[i] = sbox[s[i]]
	}
}

func ShiftRows(s []byte) {
	var t [16]byte
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			t[r+4*c] = s[r+4*((c+r)%4)]
		}
	}
	copy(s, t[:])
}

func MixColumns(s []byte) {
	for c := 0; c < 4; c++ {
		i := 4 * c
		a0, a1, a2, a3 := s[i], s[i+1], s[i+2], s[i+3]
		s[i] = xtime(a0) ^ (xtime(a1) ^ a1) ^ a2 ^ a3
		s[i+1] = a0 ^ xtime(a1) ^ (xtime(a2) ^ a2) ^ a3
		s[i+2] = a0 ^ a1 ^ xtime(a2) ^ (xtime(a3) ^ a3)
		s[i+3] = (xtime(a0) ^ a0) ^ a1 ^ a2 ^ xtime(a3)
	}
}

// Decryption round functions.

func InvSubBytes(s []byte) {
	for i := range s {
		s[i] = invSbox[s[i]]
	}
}

func InvShiftRows(s []byte) {
	var t [16]byte
	for r := 0; r < 4; r++ {
		for c := 0; c < 4; c++ {
			t[r+4*c] = s[r+4*((c-r+4)%4)]
		}
	}
	copy(s, t[:])
}

func InvMixColumns(s []byte) {
	for c := 0; c < 4; c++ {
		i := 4 * c
		a0, a1, a2, a3 := s[i], s[i+1], s[i+2], s[i+3]
		s[i] = gfMul(a0, 0x0e) ^ gfMul(a1, 0x0b) ^ gfMul(a2, 0x0d) ^ gfMul(a3, 0x09)
		s[i+1] = gfMul(a0, 0x09) ^ gfMul(a1, 0x0e) ^ gfMul(a2, 0x0b) ^ gfMul(a3, 0x0d)
		s[i+2] = gfMul(a0, 0x0d) ^ gfMul(a1, 0x09) ^ gfMul(a2, 0x0e) ^ gfMul(a3, 0x0b)
		s[i+3] = gfMul(a0, 0x0b) ^ gfMul(a1, 0x0d) ^ gfMul(a2, 0x09) ^ gfMul(a3, 0x0e)
	}
}

// Key schedule.

var rcon = [11]byte{0x00, 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36}

func KeyExpansion(key []byte) [][]byte {
	rk := make([][]byte, 11)
	rk[0] = append([]byte(nil), key...)
	for round := 1; round <= 10; round++ {
		prev := rk[round-1]
		t := [4]byte{prev[13], prev[14], prev[15], prev[12]}
		for i := range t {
			t[i] = sbox[t[i]]
		}
		t[0] ^= rcon[round]
		cur := make([]byte, 16)
		for i := 0; i < 4; i++ {
			cur[i] = prev[i] ^ t[i]
		}
		for i := 4; i < 16; i++ {
			cur[i] = prev[i] ^ cur[i-4]
		}
		rk[round] = cur
	}
	return rk
}

// Block encryption / decryption.

func EncryptBlock(pt, key []byte) []byte {
	rk := KeyExpansion(key)
	s := append([]byte(nil), pt...)

	AddRoundKey(s, rk[0])
	for round := 1; round <= 9; round++ {
		SubBytes(s)
		ShiftRows(s)
		MixColumns(s)
		AddRoundKey(s, rk[round])
	}
	SubBytes(s)
	ShiftRows(s)
	AddRoundKey(s, rk[10])
	return s
}

func DecryptBlock(ct, key []byte) []byte {
	rk := KeyExpansion(key)
	s := append([]byte(nil), ct...)

	AddRoundKey(s, rk[10])
	for round := 9; round >= 1; round-- {
		InvShiftRows(s)
		InvSubBytes(s)
		AddRoundKey(s, rk[round])
		InvMixColumns(s)
	}
	InvShiftRows(s)
	InvSubBytes(s)
	AddRoundKey(s, rk[0])
	return s
}
