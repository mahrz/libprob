#include <iostream>

#define PROB_EXPERIMENTAL
#define PROB_REDUNDANT_INFORMATION

#include "prob.hpp"

using namespace std;

RVAR_STATIC(W, 5)
RVAR_STATIC(A, 3)
RVAR_STATIC(S, 2)

RVAR(X)
RVAR(Y)
RVAR(Z)

int main(int argc, char * const argv[])
{
  prob::distribution<double, W, A, prob::given, S, W> p;
  prob::distribution<double, Z, X, prob::given, Y> q(Z(3), X(5) | Y(4));
  prob::distribution<double, W, A> r;

  W w(3), wn(1);
  A a(2);
  S s(0);

  p(w, a | s, wn) = 0.5;
  wn = w + 2;
  wn = wn - 1;
  p(w, a | s, wn) = 0.75;

  r(w, a) = 0.5;

  r.each_index([&r] (W w, A a)
  { cout << "Value here " << w << " " << a << endl;});

  p.each_index([&p] (W w, A a, prob::given, S s, W wp)
  { cout << "Value here " << w << " " << a << "|"
    << s << " " << wp << endl;});

  cout << p(W(3), A(2) | S(0), W(1)) << endl;
  cout << p(w, a | s, wn) << endl;

  cout << p << endl << endl;

  prob::distribution<double, W, A, prob::given, S, W> p1(p);
  cout << p1 << endl;

  cout << p.row_extent<1>() << endl;

  cout << r << endl << endl;

  //q.setRandom();

  q(Z(1), X(1) | Y(2)) = 0.75;
  q(Z(1), X(0) | Y(2)) = 0.75;

  q(Z(2), X(0) | Y(1)) = 3;
  q(Z(1), X(1) | Y(1)) = 2;
  q(Z(1), X(0) | Y(1)) = 4;

  q.normalize();

  cout << q << endl << endl;

  q.reshape(Z(6), X(2) | Y(3));

  cout << q << endl;

  prob::distribution<double, Z, prob::given, Y> qgrp = q.grouped_sum<0, 2>();
  cout << endl << "Marginal" << qgrp << endl;

  prob::distribution<double, X, Z, prob::given, Y> qgrp2 =
      q.grouped_sum<1, 0, 2>();

  cout << q(Z(1), X(1) | Y(2)) << endl << endl;

  cout << q.posterior_distribution(Y(2)) << endl;
  cout << q.posterior_distribution(Y(1)) << endl;

  q.each_conditional_index([&q] (Y y)
  { cout << "Condition " << y << ":" << q.posterior_distribution(y) << endl;});

  cout << q.map_copy([] (double p)
  { return 2*p;}).sum() << endl;

  prob::distribution<double, Z, X, prob::given, Y> t = q.map_copy([] (double p)
  { return sqrt(p);});

  prob::distribution<double, Y> qsum = q.sum_by_conditional();
  cout << endl << qsum << endl;

  t.map_by_conditional([] (prob::distribution<double, Z, X> post)
  { return post/post.sum();});
  cout << endl << t << endl;

  prob::distribution<double, W, prob::given, S> p_wgs;
  prob::distribution<double, W, prob::given, A> p_wga;
  prob::distribution<double, A, prob::given, S> p_ags;
  prob::distribution<double, W, prob::given, A, S> p_wgas;
  prob::distribution<double, S> p_s;
  prob::distribution<double, W> p_w;

  std::random_device rd;
  std::mt19937 gen(rd());

  prob::init::random(p_wgs, gen);
  prob::init::random(p_ags, gen);
  prob::init::random(p_wgas, gen);
  prob::init::random(p_wga, gen);
  prob::init::random(p_s, gen);
  prob::init::random(p_w, gen);

  cout << "Entropy: H(W) = " << prob::it::entropy(p_w) << endl;
  cout << "Conditional entropy: H(A|S) = "
      << prob::it::conditional_entropy(p_ags, p_s) << endl << endl;
  cout << "Mutual information: I(A;S) = "
      << prob::it::mutual_information(p_ags, p_s) << endl << endl;

  prob::distribution<double, A> p_a = uncondition(p_ags, p_s).marginalize<0>();

  cout << "Mutual information: I(A;S) = "
      << prob::it::mutual_information(p_ags, p_a, p_s) << endl << endl;

  prob::distribution<double, A> q_a;
  prob::init::random(q_a, gen);

  cout << "Div KL: D(P||Q) = " << prob::it::kl_divergence(p_a, q_a) << endl
      << endl;
  cout << "Div JS: D(P||Q) = " << prob::it::js_divergence(p_a, q_a) << endl
      << endl;

  prob::distribution<double, W, A, prob::given, S> p_wags =
      prob::join_conditionals(p_wgs, p_ags);
  std::cout << "Join Test:" << endl << p_wgs << endl << p_ags << endl << p_wags
      << endl;

  prob::distribution<double, W, S> p_ws = prob::uncondition(p_wgs, p_s);
  std::cout << "Uncondition Test:" << endl << p_wgs << endl << p_s << endl
      << p_ws << endl;

  prob::distribution<double, W, A, prob::given, S> q_wags =
      prob::partial_uncondition(p_wgas, p_ags);
  std::cout << "Partial Uncondition Test:" << endl << p_wgas << endl << p_ags
      << endl << q_wags << endl;

  prob::distribution<double, W, prob::given, S> q_wgs;
  prob::condition(p_ws, p_s, q_wgs);
  std::cout << "Condition Test:" << endl << q_wgs << endl << p_s << endl << p_ws
      << endl;

  prob::condition_conditionals(p_wags, p_ags, p_wgas);
  std::cout << "Condition Conditionals Test:" << endl << q_wags << endl << p_ags
      << endl << p_wgas << endl;

  prob::distribution<double, S, prob::given, W> p_sgw = prob::bayes(p_wgs, p_w,
      p_s);
  std::cout << "Bayes Test:" << endl << p_sgw << endl << p_wgs << endl << p_s
      << endl;

  std::cout << p_w.summary() << std::endl;

  prob::distribution<double, Z, prob::given, X, Y> pZgXY(Z(3) | X(2), Y(4));
  prob::distribution<double, X, Y> pXY(X(2), Y(4));

  prob::init::random(pZgXY, gen);
  prob::init::random(pXY, gen);

  auto pZXY = uncondition(pZgXY, pXY);
  auto pX = pZXY.marginalize<1>();
  auto pY = pZXY.marginalize<2>();
  auto pZ = pZXY.marginalize<0>();

  auto pXYZ = pZXY.marginalize<1, 2, 0>();

  prob::distribution<double, X, Y, prob::given, Z> pXYgZ(X(2), Y(4) | Z(3));
  prob::condition(pXYZ, pZ, pXYgZ);

  prob::distribution<double, Z, prob::given, X> pZgX(Z(3) | X(2));
  prob::condition(pZXY.marginalize<0, 1>(), pX, pZgX);

  prob::distribution<double, Z, prob::given, Y> pZgY(Z(3) | Y(4));
  prob::condition(pZXY.marginalize<0, 2>(), pY, pZgY);

  auto pXgZ = bayes(pZgX, pZ, pX);
  auto pYgZ = bayes(pZgY, pZ, pY);

  std::cout << prob::it::decomp::redundancy(pZgX, pZgY, pX, pY, pZ) << std::endl;
  std::cout << prob::it::decomp::minimal_information(pXgZ, pYgZ, pZ) << std::endl;
  std::cout << prob::it::decomp::icm_information(pXYgZ, pZ, pX, pY) << std::endl;

}
