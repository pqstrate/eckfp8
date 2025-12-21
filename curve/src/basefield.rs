use p3_field::extension::BinomialExtensionField;
use p3_koala_bear::KoalaBear;

/// KoalaBear degree-8 extension field
pub type BaseField = BinomialExtensionField<KoalaBear, 8>;

/// Helper function to construct a BaseField from its coefficients
#[inline]
pub fn from_coeffs(coeffs: [KoalaBear; 8]) -> BaseField {
    unsafe { core::mem::transmute(coeffs) }
}
