#import "metal_accel.h"

#if defined(__APPLE__)

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

@interface MetalAccelObj : NSObject
@property(nonatomic, strong) id<MTLDevice> device;
@property(nonatomic, strong) id<MTLCommandQueue> queue;
@property(nonatomic, strong) id<MTLComputePipelineState> pso;
@property(nonatomic, strong) id<MTLBuffer> inBuf;
@property(nonatomic, strong) id<MTLBuffer> outBuf;
@property(nonatomic, strong) id<MTLBuffer> nBuf;
@property(nonatomic, strong) id<MTLBuffer> gBuf;
@property(nonatomic, strong) id<MTLBuffer> epsBuf;
@property(nonatomic) NSUInteger capacity;
@end

@implementation MetalAccelObj
@end

struct MetalAccel {
  MetalAccelObj *obj;
};

static NSString *kKernelSrc(void) {
  return @
      "#include <metal_stdlib>\n"
      "using namespace metal;\n"
      "kernel void accel(\n"
      "  device const float4* inBodies [[buffer(0)]],\n"
      "  device float4* outAccel [[buffer(1)]],\n"
      "  constant uint &n [[buffer(2)]],\n"
      "  constant float &G [[buffer(3)]],\n"
      "  constant float &eps2 [[buffer(4)]],\n"
      "  uint id [[thread_position_in_grid]]\n"
      "){\n"
      "  if (id >= n) return;\n"
      "  float3 pi = inBodies[id].xyz;\n"
      "  float3 a = float3(0.0);\n"
      "  for (uint j = 0; j < n; j++) {\n"
      "    if (j == id) continue;\n"
      "    float3 d = inBodies[j].xyz - pi;\n"
      "    float r2 = dot(d,d) + eps2;\n"
      "    float inv = rsqrt(r2);\n"
      "    float inv3 = inv * inv * inv;\n"
      "    a += (G * inBodies[j].w) * d * inv3;\n"
      "  }\n"
      "  outAccel[id] = float4(a, 0.0);\n"
      "}\n";
}

static bool ensure_capacity(MetalAccelObj *o, NSUInteger n) {
  if (n <= o.capacity && o.inBuf && o.outBuf) return true;
  o.capacity = (n < 1024) ? 1024 : n;
  o.inBuf = [o.device newBufferWithLength:o.capacity * 4 * sizeof(float)
                                 options:MTLResourceStorageModeShared];
  o.outBuf = [o.device newBufferWithLength:o.capacity * 4 * sizeof(float)
                                  options:MTLResourceStorageModeShared];
  if (!o.inBuf || !o.outBuf) return false;
  return true;
}

MetalAccel *metal_accel_create(void) {
  @autoreleasepool {
    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev) return NULL;
    id<MTLCommandQueue> q = [dev newCommandQueue];
    if (!q) return NULL;

    NSError *err = nil;
    id<MTLLibrary> lib = [dev newLibraryWithSource:kKernelSrc() options:nil error:&err];
    if (!lib) {
      NSLog(@"Metal library compile failed: %@", err);
      return NULL;
    }
    id<MTLFunction> fn = [lib newFunctionWithName:@"accel"];
    if (!fn) return NULL;

    id<MTLComputePipelineState> pso = [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!pso) {
      NSLog(@"Metal pipeline failed: %@", err);
      return NULL;
    }

    MetalAccelObj *o = [MetalAccelObj new];
    o.device = dev;
    o.queue = q;
    o.pso = pso;
    o.capacity = 0;
    o.nBuf = [dev newBufferWithLength:sizeof(uint32_t) options:MTLResourceStorageModeShared];
    o.gBuf = [dev newBufferWithLength:sizeof(float) options:MTLResourceStorageModeShared];
    o.epsBuf = [dev newBufferWithLength:sizeof(float) options:MTLResourceStorageModeShared];
    if (!o.nBuf || !o.gBuf || !o.epsBuf) return NULL;

    MetalAccel *ma = (MetalAccel *)calloc(1, sizeof(MetalAccel));
    if (!ma) return NULL;
    ma->obj = o;
    return ma;
  }
}

void metal_accel_destroy(MetalAccel *ma) {
  if (!ma) return;
  ma->obj = nil;
  free(ma);
}

bool metal_accel_compute(MetalAccel *ma, const float *pos_mass4, size_t count, float G, float softening,
                         float *out_accel4) {
  if (!ma || !ma->obj || count == 0) return false;
  @autoreleasepool {
    MetalAccelObj *o = ma->obj;
    if (!ensure_capacity(o, (NSUInteger)count)) return false;

    memcpy([o.inBuf contents], pos_mass4, count * 4 * sizeof(float));
    *(uint32_t *)[o.nBuf contents] = (uint32_t)count;
    *(float *)[o.gBuf contents] = G;
    const float eps2 = softening * softening;
    *(float *)[o.epsBuf contents] = eps2;

    id<MTLCommandBuffer> cb = [o.queue commandBuffer];
    if (!cb) return false;
    id<MTLComputeCommandEncoder> enc = [cb computeCommandEncoder];
    if (!enc) return false;
    [enc setComputePipelineState:o.pso];
    [enc setBuffer:o.inBuf offset:0 atIndex:0];
    [enc setBuffer:o.outBuf offset:0 atIndex:1];
    [enc setBuffer:o.nBuf offset:0 atIndex:2];
    [enc setBuffer:o.gBuf offset:0 atIndex:3];
    [enc setBuffer:o.epsBuf offset:0 atIndex:4];

    const NSUInteger tg = MIN((NSUInteger)256, o.pso.maxTotalThreadsPerThreadgroup);
    MTLSize tgSize = MTLSizeMake(tg, 1, 1);
    MTLSize grid = MTLSizeMake((NSUInteger)count, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tgSize];
    [enc endEncoding];
    [cb commit];
    [cb waitUntilCompleted];

    memcpy(out_accel4, [o.outBuf contents], count * 4 * sizeof(float));
    return true;
  }
}

#else

struct MetalAccel {
  int dummy;
};

MetalAccel *metal_accel_create(void) { return NULL; }
void metal_accel_destroy(MetalAccel *ma) {
  (void)ma;
}
bool metal_accel_compute(MetalAccel *ma, const float *pos_mass4, size_t count, float G, float softening,
                         float *out_accel4) {
  (void)ma;
  (void)pos_mass4;
  (void)count;
  (void)G;
  (void)softening;
  (void)out_accel4;
  return false;
}

#endif
