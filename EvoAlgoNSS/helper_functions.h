#pragma once

template <typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
	if (!v.empty()) {
		out << '[';
		std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
		out << "\b\b]";
	}
	return out;
}

template<typename T>
auto comparator = [&](const std::vector<T>& l, const std::vector<T>& r)
{
	return f(l) < f(r);
};