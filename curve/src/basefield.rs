use p3_field::extension::BinomialExtensionField;
use p3_field::RawDataSerializable;
use p3_koala_bear::KoalaBear;

/// KoalaBear degree-8 extension field.
pub type BaseField = BinomialExtensionField<KoalaBear, 8>;

/// Helper function to construct a BaseField from its coefficients.
#[inline]
pub fn from_coeffs(coeffs: [KoalaBear; 8]) -> BaseField {
    unsafe { core::mem::transmute(coeffs) }
}

/// Serialize a BaseField element using the Plonky3 raw-data format.
pub fn to_bytes(elem: BaseField) -> [u8; <BaseField as RawDataSerializable>::NUM_BYTES] {
    let bytes: Vec<u8> = elem.into_bytes().into_iter().collect();
    bytes.try_into().expect("basefield byte length")
}

/// Serialize a BaseField element into little-endian u32 words.
pub fn to_u32s(elem: BaseField) -> [u32; <BaseField as RawDataSerializable>::NUM_BYTES / 4] {
    let words: Vec<u32> = BaseField::into_u32_stream([elem]).into_iter().collect();
    words.try_into().expect("basefield u32 length")
}
