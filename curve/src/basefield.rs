use p3_field::extension::BinomialExtensionField;
use p3_koala_bear::KoalaBear;

/// KoalaBear degree-8 extension field
pub type BaseField = BinomialExtensionField<KoalaBear, 8>;
