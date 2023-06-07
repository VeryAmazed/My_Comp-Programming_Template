//template for compare functor
struct cmp {
	bool operator()(const Class& x, const Class& y) const { return x.a < y.a; }
};