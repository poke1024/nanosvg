/*
 * TÖVE - Animated vector graphics for LÖVE.
 * https://github.com/poke1024/Tove
 *
 * Portions of this code are taken from NanoSVG.
 * Portions of this code Copyright (c) 2019, Bernhard Liebl.
 *
 * All rights reserved.
 */

void tove_deleteRasterizer(NSVGrasterizer* r) {

	if (r->stencil.data) free(r->stencil.data);
	if (r->dither.data) free(r->dither.data);
}

void tove__scanlineBit(
	NSVGrasterizer* r,
	int x,
	int y,
	int count,
	float tx,
	float ty,
	float scale,
	NSVGcachedPaint* cache,
	TOVEclip *clip) {

	unsigned char* const row = &r->bitmap[y * r->stride];
	unsigned char* const cover = r->scanline;

	const int x1 = x + count;
	for (; x < x1; x++) {
		row[x / 8] |= (cover[x] > 0 ? 1 : 0) << (x % 8);
	}
}

inline void maskClip(
	NSVGrasterizer* r,
	TOVEclip* clip,
	int xmin,
	int y,
	int count) {

	unsigned char* const cover = r->scanline;
	const int xmax = xmin + count - 1;

	for (int i = 0; i < clip->count; i++) {
		unsigned char* stencil = &r->stencil.data[
			r->stencil.size * clip->index[i] + y * r->stencil.stride];

		for (int j = xmin; j <= xmax; j++) {
			if (((stencil[j / 8] >> (j % 8)) & 1) == 0) {
				cover[j] = 0;
			}
		}
	}
}


class UnrestrictedPalette {
public:
	inline UnrestrictedPalette(
		const NSVGrasterizer *r) {
	}

	inline void operator()(int &cr, int &cg, int &cb) const {
	}
};

class RestrictedPalette {
	const TOVErasterizerQuality * const quality;

public:
	inline RestrictedPalette(
		const NSVGrasterizer *r) :
		
		quality(&r->quality) {
	}

	inline void operator()(int &cr, int &cg, int &cb) const {
		/*cr = ((cr + 16) >> 5) << 5;
		cg = ((cg + 16) >> 5) << 5;
		cb = ((cb + 16) >> 5) << 5;
		return;*/

		const int n = quality->palette.size;
		const uint8_t *p = quality->palette.colors;

		long best_d = std::numeric_limits<long>::max();
		const uint8_t *best_p = p;

		const int r0 = cr;
		const int g0 = cg;
		const int b0 = cb;

		// brute force search. fast enough for small sets.
		for (int i = 0; i < n; i++, p += 3) {
			const long dr = r0 - p[0];
			const long dg = g0 - p[1];
			const long db = b0 - p[2];

			const long d = dr * dr + dg * dg + db * db;
			if (d < best_d) {
				best_d = d;
				best_p = p;
			}
		}

		cr = best_p[0];
		cg = best_p[1];
		cb = best_p[2];
	}
};

typedef UnrestrictedPalette AnyPalette;


class NoDithering {
public:
	inline NoDithering(
		NSVGrasterizer* r,
		NSVGcachedPaint *cache,
		int x,
		int y,
		int count) {
	}

	template<typename Palette>
	inline void operator()(
		int r0, int g0, int b0, int a0,
		int &r, int &g, int &b, int& a,
		const int x,
		const Palette &palette) {

		r = r0;
		g = g0;
		b = b0;
		a = a0;

		palette(r, g, b);
	}
};

class SierraDithering {
public:
	typedef float dither_error_t;
	
	static constexpr int diffusion_matrix_half_w = 2;
	static constexpr int diffusion_matrix_height = 3;
	static constexpr int dither_components = 4;

protected:
	dither_error_t *diffusion[diffusion_matrix_height];

	inline void rotate(
		const NSVGrasterizer* r,
		NSVGcachedPaint* cache,
		const int x,
		const int y,
		const int count) {

		dither_error_t **rows = diffusion;
		dither_error_t * const data = static_cast<dither_error_t*>(r->dither.data);

		const unsigned int stride = r->dither.stride;

		int roll = y - cache->tove.ditherY;
		if (roll < 0) { // should never happen.
			roll = diffusion_matrix_height;
		} else if (roll > diffusion_matrix_height) {
			roll = diffusion_matrix_height;
		}

		cache->tove.ditherY = y;

		// rotate rows.
		for (int i = 0; i < diffusion_matrix_height; i++) {
			rows[i] = &data[((y + i) % diffusion_matrix_height) * stride +
				diffusion_matrix_half_w * dither_components];
		}

		const int span0 = (x - diffusion_matrix_half_w) * dither_components;
		const int span1 = (x + count + diffusion_matrix_half_w) * dither_components;

		// new rows.
		for (int i = 0; i < roll; i++) {
			dither_error_t * const row = rows[diffusion_matrix_height - 1 - i];
			for (int j = span0; j < span1; j++) {
				row[j] = 0.0f;
			}
		}

		// old rows that already have content.
		for (int i = 0; i < diffusion_matrix_height - roll; i++) {
			dither_error_t * const row = rows[i];

			const int previousLeftSpan = cache->tove.ditherSpan[0];
			for (int j = span0; j < previousLeftSpan; j++) {
				row[j] = 0.0f;
			}
			
			for (int j = cache->tove.ditherSpan[1]; j < span1; j++) {
				row[j] = 0.0f;
			}
		}

#if 0 // visual results of this are quite mixed.
		if (roll < diffusion_matrix_height) {
			// distribute errors from previous rows from the left.

			dither_error_t * const r0 = rows[0];
			dither_error_t * const r1 = rows[1];
			dither_error_t * const r2 = rows[2];

			for (int k = cache->tove.ditherSpan[0] / dither_components; k < x; k++) {
				const float f_cr = r0[k * dither_components + 0];
				const float f_cg = r0[k * dither_components + 1];
				const float f_cb = r0[k * dither_components + 2];
				const float f_ca = r0[k * dither_components + 3];

				const float errors[] = {f_cr, f_cg, f_cb, f_ca};

				if (errors[0] == 0.0f &&
					errors[1] == 0.0f &&
					errors[2] == 0.0f &&
					errors[3] == 0.0f) {
					continue;
				}

				// distribute errors using sierra dithering.
				#pragma clang loop vectorize(enable)
				for (int j = 0; j < 4; j++) {
					const int offset = k * dither_components + j;
					const float error = errors[j];

					r0[offset + 1 * dither_components] += error * (5.0f / 32.0f);
					r0[offset + 2 * dither_components] += error * (3.0f / 32.0f);

					r1[offset - 2 * dither_components] += error * (2.0f / 32.0f);
					r1[offset - 1 * dither_components] += error * (4.0f / 32.0f);
					r1[offset + 0 * dither_components] += error * (5.0f / 32.0f);
					r1[offset + 1 * dither_components] += error * (4.0f / 32.0f);
					r1[offset + 2 * dither_components] += error * (2.0f / 32.0f);

					r2[offset - 1 * dither_components] += error * (2.0f / 32.0f);
					r2[offset + 0 * dither_components] += error * (3.0f / 32.0f);
					r2[offset + 1 * dither_components] += error * (2.0f / 32.0f);
				}
			}
		}
#endif

		cache->tove.ditherSpan[0] = span0;
		cache->tove.ditherSpan[1] = span1;
	}

public:
	static inline bool enabled(const NSVGrasterizer* r) {
		return r->quality.flags > 0;
	}

	static inline bool allocate(NSVGrasterizer* r, int w) {
		if (enabled(r)) {
			TOVEdither &dither = r->dither;
			dither.stride = (w + 2 * diffusion_matrix_half_w) * dither_components;
			dither.data = (dither_error_t*)realloc(dither.data,
				dither.stride * diffusion_matrix_height * sizeof(dither_error_t));
			return dither.data != nullptr;
		} else {
			return true;
		}
	}

	inline SierraDithering(
		NSVGrasterizer* r,
		NSVGcachedPaint *cache,
		int x,
		int y,
		int count) {

		rotate(r, cache, x, y, count);
	}

	template<typename Palette>
	inline void operator()(
		int r0, int g0, int b0, int a0,
		int &r, int &g, int &b, int& a,
		const int x,
		const Palette &palette) {

		const uint32_t color = (*this)(
			r0,
			g0,
			b0,
			a0,
			x,
			palette);
		
		r = color & 0xff;
		g = (color >> 8) & 0xff;
		b = (color >> 16) & 0xff;
		a = (color >> 24) & 0xff;
	}

	template<typename Palette>
	inline uint32_t operator()(
		float r, float g, float b, float a,
		const int x,
		const Palette &palette) {
		
		dither_error_t * const r0 = diffusion[0];
		dither_error_t * const r1 = diffusion[1];
		dither_error_t * const r2 = diffusion[2];

		const float f_cr = nsvg__clampf(r + r0[x * dither_components + 0], 0.0f, 255.0f);
		const float f_cg = nsvg__clampf(g + r0[x * dither_components + 1], 0.0f, 255.0f);
		const float f_cb = nsvg__clampf(b + r0[x * dither_components + 2], 0.0f, 255.0f);
		const float f_ca = nsvg__clampf(a + r0[x * dither_components + 3], 0.0f, 255.0f);

		int cr = f_cr + 0.5f;
		int cg = f_cg + 0.5f;
		int cb = f_cb + 0.5f;
		const int ca = f_ca + 0.5f;

		palette(cr, cg, cb);

		const float errors[] = {f_cr - cr, f_cg - cg, f_cb - cb, f_ca - ca};

		// distribute errors using sierra dithering.
		#pragma clang loop vectorize(enable)
		for (int j = 0; j < 4; j++) {
			const int offset = x * dither_components + j;
			const float error = errors[j];

			r0[offset + 1 * dither_components] += error * (5.0f / 32.0f);
			r0[offset + 2 * dither_components] += error * (3.0f / 32.0f);

			r1[offset - 2 * dither_components] += error * (2.0f / 32.0f);
			r1[offset - 1 * dither_components] += error * (4.0f / 32.0f);
			r1[offset + 0 * dither_components] += error * (5.0f / 32.0f);
			r1[offset + 1 * dither_components] += error * (4.0f / 32.0f);
			r1[offset + 2 * dither_components] += error * (2.0f / 32.0f);

			r2[offset - 1 * dither_components] += error * (2.0f / 32.0f);
			r2[offset + 0 * dither_components] += error * (3.0f / 32.0f);
			r2[offset + 1 * dither_components] += error * (2.0f / 32.0f);
		}

		/*cr = 128 + 0.5 * nsvg__clampf(errors[0],-255, 255);
		cg = 128 + 0.5 * nsvg__clampf(errors[0],-255, 255);
		cb = 128 + 0.5 * nsvg__clampf(errors[0],-255, 255);*/

		return cr | (cg << 8) | (cb << 16) | (ca << 24);
	}
};


class LinearGradient {
	const float* const t;

public:
	inline LinearGradient(const NSVGcachedPaint* cache) : t(cache->xform) {
	}

	inline float operator()(float fx, float fy) const {
		return fx*t[1] + fy*t[3] + t[5];
	}
};

class RadialGradient {
	const float* const t;

public:
	inline RadialGradient(const NSVGcachedPaint* cache) : t(cache->xform) {
	}

	inline float operator()(float fx, float fy) const {
		const float gx = fx*t[0] + fy*t[2] + t[4];
		const float gy = fx*t[1] + fy*t[3] + t[5];
		return sqrtf(gx*gx + gy*gy);
	}
};


class FastGradientColors {
	const NSVGcachedPaint * const cache;

public:
	inline FastGradientColors(
		NSVGrasterizer* r,
		NSVGcachedPaint *cache,
		int x,
		int y,
		int count) : cache(cache) {
	}

	inline bool good() {
		return true;
	}

	inline uint32_t operator()(int x, float gy) const {
		return cache->colors[(int)nsvg__clampf(gy*255.0f, 0, 255.0f)];
	}
};

template<typename Dithering, typename Palette>
class BestGradientColors {

	const NSVGgradient * const gradient;
	const NSVGgradientStop *stop;
	const NSVGgradientStop *const stopN;

	Dithering dithering;
	const Palette palette;

public:
	inline BestGradientColors(
		NSVGrasterizer* r,
		NSVGcachedPaint *cache,
		int x,
		int y,
		int count) :

		dithering(r, cache, x, y, count),
		palette(r),

		gradient(cache->tove.paint->gradient),
		stop(gradient->stops),
		stopN(gradient->stops + gradient->nstops) {

	}

	static void init(
		NSVGcachedPaint* cache,
		NSVGpaint* paint,
		float opacity) {

		cache->tove.paint = paint;
		cache->tove.ditherY = -Dithering::diffusion_matrix_height;

		NSVGgradientStop *stop = paint->gradient->stops;
		const int n = paint->gradient->nstops;

		for (int i = 0; i < n; i++) {
			stop->tove.color = nsvg__applyOpacity(stop->color, opacity);
			stop->tove.offset = nsvg__clampf(stop->offset, 0.0f, 1.0f);
			stop++;
		}
	}

	inline bool good() const {
		return stop + 1 < stopN;
	}

	inline uint32_t operator()(const int x, const float gy) {
		while (gy >= (stop + 1)->tove.offset && (stop + 2) < stopN) {
			stop++;
		}
		while (gy < stop->tove.offset && stop > gradient->stops) {
			stop--;
		}

		const unsigned int c0 = stop->tove.color;
		const int cr0 = (c0) & 0xff;
		const int cg0 = (c0 >> 8) & 0xff;
		const int cb0 = (c0 >> 16) & 0xff;
		const int ca0 = (c0 >> 24) & 0xff;

		const unsigned int c1 = (stop + 1)->tove.color;
		const int cr1 = (c1) & 0xff;
		const int cg1 = (c1 >> 8) & 0xff;
		const int cb1 = (c1 >> 16) & 0xff;
		const int ca1 = (c1 >> 24) & 0xff;

		const float offset0 = stop->tove.offset;
		const float range = (stop + 1)->tove.offset - offset0;
		const float t = nsvg__clampf((gy - offset0) / range, 0.0f, 1.0f);
		const float s = 1.0f - t;

		return dithering(
			cr0 * s + cr1 * t,
			cg0 * s + cg1 * t,
			cb0 * s + cb1 * t,
			ca0 * s + ca1 * t,
			x,
			palette);
	}
};


template<typename Gradient, typename Colors>
void drawGradientScanline(
	NSVGrasterizer* r,
	int x,
	int y,
	int count,
	float tx,
	float ty,
	float scale,
	NSVGcachedPaint* cache,
	TOVEclip* clip) {

	using tove::nsvg__div255;

	unsigned char* dst0 = &r->bitmap[y * r->stride] + x*4;
	unsigned char* cover0 = &r->scanline[x];
	maskClip(r, clip, x, y, count);

	// TODO: spread modes.
	// TODO: plenty of opportunities to optimize.
	float fx, fy, dx, gy;
	int i, cr, cg, cb, ca;
	unsigned int c;
	const Gradient gradient(cache);
	Colors colors(r, cache, x, y, count);

	fx = ((float)x - tx) / scale;
	fy = ((float)y - ty) / scale;
	dx = 1.0f / scale;

	// use serpentine scanning for better dithering.
	int i0, id;
	if (y & 1) {
		i0 = 0;
		id = 1;
	} else {
		fx += (count - 1) * dx;
		dx = -dx;
		i0 = count - 1;
		id = -1;
	}

	for (int k = 0; k < count; k++) {
		const int i = i0 + k * id;

		int r,g,b,a,ia;
		c = colors(x + i, gradient(fx, fy));
		cr = (c) & 0xff;
		cg = (c >> 8) & 0xff;
		cb = (c >> 16) & 0xff;
		ca = (c >> 24) & 0xff;

		unsigned char *cover = cover0 + i;
		a = nsvg__div255((int)cover[0] * ca);
		ia = 255 - a;

		// Premultiply
		r = nsvg__div255(cr * a);
		g = nsvg__div255(cg * a);
		b = nsvg__div255(cb * a);

		// Blend over
		unsigned char *dst = dst0 + 4 * i;
		r += nsvg__div255(ia * (int)dst[0]);
		g += nsvg__div255(ia * (int)dst[1]);
		b += nsvg__div255(ia * (int)dst[2]);
		a += nsvg__div255(ia * (int)dst[3]);

		dst[0] = (unsigned char)r;
		dst[1] = (unsigned char)g;
		dst[2] = (unsigned char)b;
		dst[3] = (unsigned char)a;

		fx += dx;
	}
}

template<typename Dithering, typename Palette>
inline void drawColorScanline(
	NSVGrasterizer* r,
	int xmin,
	int y,
	int count,
	float tx,
	float ty,
	float scale,
	NSVGcachedPaint* cache,
	TOVEclip* clip) {

	unsigned char* dst = &r->bitmap[y * r->stride] + xmin*4;
	unsigned char* cover = &r->scanline[xmin];
	maskClip(r, clip, xmin, y, count);

	const int cr0 = cache->colors[0] & 0xff;
	const int cg0 = (cache->colors[0] >> 8) & 0xff;
	const int cb0 = (cache->colors[0] >> 16) & 0xff;
	const int ca0 = (cache->colors[0] >> 24) & 0xff;

	Dithering dithering(r, cache, xmin, y, count);
	const Palette palette(r);

	for (int i = 0; i < count; i++) {
		int cr, cg, cb, ca;
		dithering(
			cr0, cg0, cb0, ca0,
			cr, cg, cb, ca,
			xmin + i,
			palette);

		int r,g,b;
		int a = nsvg__div255((int)cover[0] * ca);
		int ia = 255 - a;
		// Premultiply
		r = nsvg__div255(cr * a);
		g = nsvg__div255(cg * a);
		b = nsvg__div255(cb * a);

		// Blend over
		r += nsvg__div255(ia * (int)dst[0]);
		g += nsvg__div255(ia * (int)dst[1]);
		b += nsvg__div255(ia * (int)dst[2]);
		a += nsvg__div255(ia * (int)dst[3]);

		dst[0] = (unsigned char)r;
		dst[1] = (unsigned char)g;
		dst[2] = (unsigned char)b;
		dst[3] = (unsigned char)a;

		cover++;
		dst += 4;
	}
}

TOVEscanlineFunction tove__initPaint(
	NSVGcachedPaint* cache,
	const NSVGrasterizer* r,
	NSVGpaint* paint,
	float opacity,
	bool &initCacheColors) {

	if (r && SierraDithering::enabled(r)) {
		switch (cache->type) {
			case NSVG_PAINT_COLOR:
				initCacheColors = true;
				if (r->quality.palette.colors) {
					return drawColorScanline<SierraDithering, RestrictedPalette>;
				} else {
					return drawColorScanline<NoDithering, UnrestrictedPalette>;
				}
				break;
			case NSVG_PAINT_LINEAR_GRADIENT:
				BestGradientColors<SierraDithering, AnyPalette>::init(cache, paint, opacity);
				initCacheColors = false;
				if (r->quality.palette.colors) {
					return drawGradientScanline<
						LinearGradient,
						BestGradientColors<SierraDithering, RestrictedPalette>>;
				} else {
					return drawGradientScanline<
						LinearGradient,
						BestGradientColors<SierraDithering, UnrestrictedPalette>>;
				}
				break;
			case NSVG_PAINT_RADIAL_GRADIENT:
				BestGradientColors<SierraDithering, AnyPalette>::init(cache, paint, opacity);
				initCacheColors = false;
				if (r->quality.palette.colors) {
					return drawGradientScanline<
						RadialGradient,
						BestGradientColors<SierraDithering, RestrictedPalette>>;
				} else {
					return drawGradientScanline<
						RadialGradient,
						BestGradientColors<SierraDithering, UnrestrictedPalette>>;
				}
				break;
		}
	}

	switch (cache->type) {
		case NSVG_PAINT_LINEAR_GRADIENT:
			initCacheColors = true;
			return drawGradientScanline<LinearGradient, FastGradientColors>;
		case NSVG_PAINT_RADIAL_GRADIENT:
			initCacheColors = true;
			return drawGradientScanline<RadialGradient, FastGradientColors>;
		default:
			initCacheColors = false;
			return drawColorScanline<NoDithering, UnrestrictedPalette>;
	}
}

bool tove__rasterize(
	NSVGrasterizer* r,
    NSVGimage* image,
    int w,
    int h,
	float tx,
    float ty,
    float scale)
{
	if (!SierraDithering::allocate(r, w)) {
		return false;
	}

	TOVEclipPath* clipPath;
	int clipPathCount = 0;

	clipPath = image->clip.instances;
	if (clipPath == NULL) {
		return true;
	}

	while (clipPath != NULL) {
		clipPathCount++;
		clipPath = clipPath->next;
	}

	r->stencil.stride = w / 8 + (w % 8 != 0 ? 1 : 0);
	r->stencil.size = h * r->stencil.stride;
	r->stencil.data = (unsigned char*)realloc(
		r->stencil.data, r->stencil.size * clipPathCount);
	if (r->stencil.data == NULL) {
		return false;
	}
	memset(r->stencil.data, 0, r->stencil.size * clipPathCount);

	clipPath = image->clip.instances;
	while (clipPath != NULL) {
		nsvg__rasterizeShapes(r, clipPath->shapes, tx, ty, scale,
			&r->stencil.data[r->stencil.size * clipPath->index],
			w, h, r->stencil.stride, tove__scanlineBit);
		clipPath = clipPath->next;
	}

	return true;
}
