using Unitful

↦(value::Unitful.Units, unit::Unitful.Units) = Unitful.convfact(unit, value)
↦(value::Unitful.Quantity, unit::Unitful.Units) = uconvert(unit, value).val

# units
const m = u"m" ↦ u"m"
const u = u"u" ↦ u"kg"
const C = u"C" ↦ u"C"
const V = u"V" ↦ u"V"
const MV = u"MV" ↦ u"V"
const eV = u"eV" ↦ u"eV"
const kg = u"kg" ↦ u"kg"
const kmps = u"km/s" ↦ u"m/s"

# constants
const ɛ0 = u"ɛ0" ↦ u"F/m"
const e  = 1.602_176_6208e-19C  # elementary charge
const me = 9.109_383_7015e-31kg # mass of electron
