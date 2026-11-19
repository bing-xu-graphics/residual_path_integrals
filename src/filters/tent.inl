Vector2 sample_op::operator()(const Tent &filter) const {

    // for the positive part, we map [0.5, 1] to [0, 1] again
    // the integral is
    // (1/(h/2)) * (x - x^2/(2h))
    // and inverting (-x^2 + 2hx- v h^2 = 0) we get
    // x = (-2h + sqrt(4h^2 - 4 v h^2))) / (-2)
    //   = h - h sqrt(1 - v)
    // when v = 0, x = 0, v = 1, x = h

    Real h = filter.width / 2;
    Real x = rnd_param[0] < 0.5 ? 
        h * (sqrt(2 * rnd_param[0]) - 1) :
        h * (1 - sqrt(1 - 2 * (rnd_param[0] - Real(0.5))));
    Real y = rnd_param[1] < 0.5 ? 
        h * (sqrt(2 * rnd_param[1]) - 1) :
        h * (1 - sqrt(1 - 2 * (rnd_param[1] - Real(0.5))));
    return Vector2{x, y};
}
