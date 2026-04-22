use ff::PrimeField;
use merlin::Transcript;
use pasta_curves::group::{Curve, GroupEncoding};
use pasta_curves::pallas::{Affine as PallasAffine, Point as PallasPoint, Scalar as PallasScalar};
use crate::traits::ByteSerializable;

pub fn append_to_transcript<C: ByteSerializable>(t: &mut Transcript, label: &'static [u8], elem: &C) {
    t.append_message(label, &elem.serialize_to_bytes());
}

pub fn get_challenge<C: ByteSerializable>(t: &mut Transcript, label: &'static [u8]) -> C {
    let mut buf = [0u8; 64];
    t.challenge_bytes(label, &mut buf);
    C::deserialize_from_bytes(&buf).unwrap()
}

pub fn get_challenge_u64(transcript: &mut Transcript, label: &'static [u8]) -> u64 {
    let mut buf = [0u8; 8];
    transcript.challenge_bytes(label, &mut buf);
    u64::from_le_bytes(buf)
}

/// Append a scalar to the transcript
pub fn transcript_append_scalar(t: &mut Transcript, label: &'static [u8], s: &PallasScalar) {
    t.append_message(label, &s.serialize_to_bytes());
}

/// Append a group element to the transcript
pub fn transcript_append_point(t: &mut Transcript, label: &'static [u8], p: &PallasPoint) {
    let _a = p.to_affine();
    t.append_message(label, &p.serialize_to_bytes());
}

/// Derive a challenge scalar from transcript
pub fn transcript_challenge_scalar(t: &mut Transcript, label: &'static [u8]) -> PallasScalar {
    let mut buf = [0u8; 64];
    t.challenge_bytes(label, &mut buf);
    // Take first 32 bytes and interpret as 4 u64 limbs (little-endian)
    let mut limbs = [0u64; 4];
    for i in 0..4 {
        limbs[i] = u64::from_le_bytes(buf[i * 8..(i + 1) * 8].try_into().unwrap());
    }
    PallasScalar::from_raw(limbs)
}
