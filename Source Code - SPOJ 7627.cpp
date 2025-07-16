#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <algorithm>
using namespace std;
 
const double eps = 1e-08;
inline int fi (const double& a)
{
	return a > eps ? 1 : (a >= -eps ? 0 : -1);
}
const double pi = 3.14159265358979;
const double hpi = pi * 0.5;
inline double sta (const double& a)
{
	if (fi(a - pi) > 0) return a - 2 * pi;
	else if (fi(a + pi) <= 0) return a + 2 * pi;
	else return a;
}
 
double r;
typedef pair<double, double> pdd;
#define X first
#define Y second
inline double len (const pdd& a, const pdd& b)
{
	return sqrt((a.X - b.X) * (a.X - b.X) + (a.Y - b.Y) * (a.Y - b.Y));
}
struct rec
{
	double x0, y0, x1, y1;
	rec (void) {}
	rec (double _x0, double _y0, double _x1, double _y1) :
		x0(_x0), y0(_y0), x1(_x1), y1(_y1) {}
};
 
void outer_tan (const pdd& a, const pdd& b, double ra[], double rb[])
{
	double sa = atan2(b.Y - a.Y, b.X - a.X), sb = sta(sa + pi);
	ra[0] = sta(sa + hpi), rb[0] = sta(sb - hpi);
	ra[1] = sta(sa - hpi), rb[1] = sta(sb + hpi);
}
bool inner_tan (const pdd& a, const pdd& b, double ra[], double rb[])
{
	double l = len(a, b);
	if (fi(l - r * 2) <= 0) return false;
	double sa = atan2(b.Y - a.Y, b.X - a.X), sb = sta(sa + pi);
	double t = acos(r * 2 / l);
	ra[0] = sta(sa + t), rb[0] = sta(sb + t);
	ra[1] = sta(sa - t), rb[1] = sta(sb - t);
	return true;
}
void hori_blk (const pdd& c, double y, double x0, double x1, vector<double>& b)
{
	double l = y - c.Y;
	if (fi(l - r) < 0 && fi(l + r) > 0)
	{
		double t = acos(l / r), u0 = c.X + r * cos(hpi - t), u1 = c.X + r * cos(hpi + t);
		if (fi(u0 - x0) > 0 && fi(u0 - x1) < 0) b.push_back(sta(hpi - t));
		if (fi(u1 - x0) > 0 && fi(u1 - x1) < 0) b.push_back(sta(hpi + t));
	}
}
void verti_blk (const pdd& c, double x, double y0, double y1, vector<double>& b)
{
	double l = x - c.X;
	if (fi(l - r) < 0 && fi(l + r) > 0)
	{
		double t = acos(l / r), u0 = c.Y + r * sin(-t), u1 = c.Y + r * sin(t);
		if (fi(u0 - y0) > 0 && fi(u0 - y1) < 0) b.push_back(sta(-t));
		if (fi(u1 - y0) > 0 && fi(u1 - y1) < 0) b.push_back(sta(t));
	}
}
void rec_blk (const rec& tr, const pdd& c, vector<double>& b)
{
	hori_blk(c, tr.y0, tr.x0, tr.x1, b), hori_blk(c, tr.y1, tr.x0, tr.x1, b);
	verti_blk(c, tr.x0, tr.y0, tr.y1, b), verti_blk(c, tr.x1, tr.y0, tr.y1, b);
}
void cir_blk (const pdd& tr, const pdd& c, vector<double>& b)
{
	double l = len(tr, c);
	if (fi(l - 2 * r) < 0)
	{
		double s = atan2(tr.Y - c.Y, tr.X - c.X), t = acos(l / (2 * r));
		b.push_back(sta(s + t)), b.push_back(sta(s - t));
	}
}
inline double cross (double x0, double y0, double x1, double y1)
{
	return x0 * y1 - x1 * y0;
}
inline double dot (double x0, double y0, double x1, double y1)
{
	return x0 * x1 + y0 * y1;
}
inline int fic (double x0, double y0, double x1, double y1)
{
	return fi(cross(x0, y0, x1, y1));
}
inline int fid (double x0, double y0, double x1, double y1)
{
	return fi(dot(x0, y0, x1, y1));
}
inline bool seg_isct (const pdd& a0, const pdd& b0, const pdd& a1, const pdd& b1)
{
	return (fic(a1.X - a0.X, a1.Y - a0.Y, a1.X - b0.X, a1.Y - b0.Y) *
			fic(b1.X - b0.X, b1.Y - b0.Y, b1.X - a0.X, b1.Y - a0.Y) == 1 && 
			fic(a0.X - a1.X, a0.Y - a1.Y, a0.X - b1.X, a0.Y - b1.Y) *
			fic(b0.X - b1.X, b0.Y - b1.Y, b0.X - a1.X, b0.Y - a1.Y) == 1);
}
inline bool rec_in (const rec& tr, const pdd& a)
{
	return (fi(tr.x0 - a.X) < 0 && fi(a.X - tr.x1) < 0 &&
			fi(tr.y0 - a.Y) < 0 && fi(a.Y - tr.y1) < 0);
}
bool rec_isct (const rec& tr, const pdd& a, const pdd& b)
{
	if (seg_isct(a, b, pdd(tr.x0, tr.y0), pdd(tr.x0, tr.y1)) ||
		seg_isct(a, b, pdd(tr.x1, tr.y0), pdd(tr.x1, tr.y1)) ||
		seg_isct(a, b, pdd(tr.x0, tr.y0), pdd(tr.x1, tr.y0)) ||
		seg_isct(a, b, pdd(tr.x0, tr.y1), pdd(tr.x1, tr.y1))) return true;
	else if (rec_in(tr, pdd((a.X + b.X) * 0.5, (a.Y + b.Y) * 0.5))) return true;
	else return false;
}
bool cir_isct (const pdd& c, const pdd& a, const pdd& b)
{
	double la = len(c, a), lb = len(c, b);
	if (fi(la - r) < 0 || fi(lb - r) < 0) return true;
	else if (fid(b.X - c.X, b.Y - c.Y, b.X - a.X, b.Y - a.Y) > 0 &&
			 fid(a.X - c.X, a.Y - c.Y, a.X - b.X, a.Y - b.Y) > 0)
	{
		double l = len(a, b), s2 = fabs(cross(a.X - c.X, a.Y - c.Y, b.X - c.X, b.Y - c.Y));
		double h = s2 / l;
		return fi(h - r) < 0;
	}
	else return false;
}
 
const int N = 262144;
const int M = 500010;
int head[N], to[M], nxt[M]; double ll[N];
int hmr, gmr;
inline int new_node (void)
{
	head[hmr] = -1;
	return hmr++;
}
inline void link (int a, int b, double l)
{
	int p = gmr++;
	to[p] = b, nxt[p] = head[a], ll[p] = l;
	head[a] = p;
}
double dist[N], td[N];
int cap, val[N * 2];
inline void make_cap (void)
{
	cap = 2; while (cap < hmr) cap <<= 1;
	for (int i = 0; i < cap; i++) val[i | cap] = i;
}
inline int min_vertex (int a, int b)
{
	if (fi(dist[a] - dist[b]) < 0) return a;
	else return b;
}
inline void modify (int a, double v)
{
	dist[a] = v, a |= cap;
	while (a != 1)
	{
		a >>= 1;
		val[a] = min_vertex(val[a << 1], val[(a << 1) | 1]);
	}
}
const double mi = 1e+100, mie = 1e+80;
double dijkstra (int source, int ter)
{
	make_cap();
	for (int i = 0; i < cap; i++) dist[i] = td[i] = mi;
	td[source] = 0, modify(source, 0);
	while (1)
	{
		if (dist[val[1]] > mie) break;
		int c = val[1]; modify(c, mi);
		for (int p = head[c]; p != -1; p = nxt[p])
		{
			int t = to[p]; double l = ll[p];
			if (fi(td[c] + l - td[t]) < 0)
			{
				td[t] = td[c] + l;
				modify(t, td[t]);
			}
		}
	}
	return td[ter];
}
 
pdd cir[160]; rec rect[80]; int rm, cm;
typedef pair<double, int> pdi;
vector<double> block[160];
vector<pdi> key[160];
bool valid_seg (const pdd& a, const pdd& b)
{
	for (int i = 0; i < rm; i++) if (rec_isct(rect[i], a, b)) return false;
	for (int i = 0; i < cm; i++) if (cir_isct(cir[i], a, b)) return false;
	return true;
}
bool valid_pt (const pdd& p)
{
    for (int i = 0; i < rm; i++) if (rec_in(rect[i], p)) return false;
    for (int i = 0; i < cm; i++) if (fi(len(p, cir[i]) - r) < 0) return false;
    return true;
}
inline void complex_link (int i, int j, double ai[], double aj[])
{
	for (int k = 0; k < 2; k++)
	{
		pdd ui(cir[i].X + cos(ai[k]) * r, cir[i].Y + sin(ai[k]) * r),
			uj(cir[j].X + cos(aj[k]) * r, cir[j].Y + sin(aj[k]) * r);
		if (valid_seg(ui, uj))
		{
			double l = len(ui, uj);
			int pi = new_node(), pj = new_node();
			key[i].push_back(pdi(ai[k], pi)), key[j].push_back(pdi(aj[k], pj));
			link(pi, pj, l), link(pj, pi, l);
		}
	} 
}
int main ()
{
	int n; scanf("%lf %d", &r, &n);
	double xt[2], yt[2], ans; scanf("%lf %lf %lf %lf", &xt[0], &yt[0], &xt[1], &yt[1]);
	for (int i = 0; i < n; i++)
	{
		double x0, y0, x1, y1; scanf("%lf %lf %lf %lf", &x0, &y0, &x1, &y1);
		rect[i << 1] = rec(x0 - r, y0, x1 + r, y1);
		rect[(i << 1) | 1] = rec(x0, y0 - r, x1, y1 + r);
		cir[i << 2] = pdd(x0, y0), cir[(i << 2) | 1] = pdd(x0, y1);
		cir[(i << 2) | 2] = pdd(x1, y1), cir[(i << 2) | 3] = pdd(x1, y0);
	}
	cm = (n << 2), rm = (n << 1);
	if (valid_seg(pdd(xt[0], yt[0]), pdd(xt[1], yt[1])))
	{
		ans = len(pdd(xt[0], yt[0]), pdd(xt[1], yt[1]));
	}
	else
	{
		new_node(), new_node();
		for (int i = 0; i < cm; i++)
		{
			for (int j = i + 1; j < cm; j++)
			{
				double ai[2], aj[2];
				outer_tan(cir[i], cir[j], ai, aj);
				complex_link(i, j, ai, aj);
				if (inner_tan(cir[i], cir[j], ai, aj))
				{
					complex_link(i, j, ai, aj);
                }
                else if (fi(len(cir[i], cir[j]) - r * 2) == 0)
                {
                    double si = atan2(cir[j].Y - cir[i].Y, cir[j].X - cir[i].X), sj = sta(si - pi);
                    if (valid_pt(pdd(cir[i].X + r * cos(si), cir[i].Y + r * sin(si))))
                    {
                        int pi = new_node(), pj = new_node();
                        key[i].push_back(pdi(si, pi)), key[j].push_back(pdi(sj, pj));
                        link(pi, pj, 0), link(pj, pi, 0);
                    }
                }
			}
		}
		for (int dt = 0; dt < 2; dt++)
		{
			double ai[2], av[2]; pdd mt(xt[dt], yt[dt]);
			for (int i = 0; i < cm; i++)
			{
				if (fi(len(mt, cir[i]) - r) == 0)
				{
					int p = new_node();
					double s = atan2(yt[dt] - cir[i].Y, xt[dt] - cir[i].X);
					key[i].push_back(pdi(s, p));
					link(dt, p, 0), link(p, dt, 0);
				}
				else
				{
					pdd v_cir(xt[dt] * 2 - cir[i].X, yt[dt] * 2 - cir[i].Y);
					bool tt = inner_tan(cir[i], v_cir, ai, av);
					for (int k = 0; k < 2; k++)
					{
						pdd u(cir[i].X + cos(ai[k]) * r, cir[i].Y + sin(ai[k]) * r);
						if (valid_seg(mt, u))
						{
							double l = len(mt, u);
							int p = new_node();
							key[i].push_back(pdi(ai[k], p));
							link(dt, p, l), link(p, dt, l);
						}
					}
				}
			}
		}
		for (int i = 0; i < cm; i++)
		{
			for (int j = 0; j < cm; j++)
			{
				if (j == i) continue;
				cir_blk(cir[j], cir[i], block[i]);
			}
			for (int j = 0; j < rm; j++)
			{
				rec_blk(rect[j], cir[i], block[i]);
			}
			sort(key[i].begin(), key[i].end());
			sort(block[i].begin(), block[i].end());
			vector<pdi>::iterator pt = key[i].begin(), fore = key[i].end();
            vector<double>::iterator bt = block[i].begin();
			while (pt != key[i].end())
			{
                if (fore == key[i].end())
                {
                    while (bt != block[i].end() && fi(*bt - pt->X) <= 0) ++bt;
                    fore = pt++;
                }
                else
                {
                    if (bt == block[i].end() || fi(pt->first - *bt) <= 0)
                    {
                        double da = pt->X - fore->X;
                        link(fore->Y, pt->Y, da * r), link(pt->Y, fore->Y, da * r);
                        fore = pt++;
                    }
                    else fore = key[i].end();
                }
            }
            int ks = key[i].size() - 1, bs = block[i].size() - 1;
            if (ks > 0 && (bs < 0 || (fi(key[i][0].X - block[i][0]) <= 0 && fi(key[i][ks].X - block[i][bs]) >= 0)))
            {
                double da = 2 * pi - (key[i][ks].X - key[i][0].X);
                link(key[i][0].Y, key[i][ks].Y, da * r), link(key[i][ks].Y, key[i][0].Y, da * r);
            }
		}
		ans = dijkstra(0, 1);
	}
	if (ans > mie) printf("no solution\n");
	else printf("%.6f\n", ans);
	return 0;
}
