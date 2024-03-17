# Most Important Step
Add shader in pubspec under shader tag , not assets 


# 1
* Reference [Link](https://www.shadertoy.com/view/mtyGWy)
```glsl
/* This animation is the material of my first youtube tutorial about creative 
   coding, which is a video in which I try to introduce programmers to GLSL 
   and to the wonderful world of shaders, while also trying to share my recent 
   passion for this community.
                                       Video URL: https://youtu.be/f4s1h2YETNY
*/

//https://iquilezles.org/articles/palettes/
vec3 palette( float t ) {
    vec3 a = vec3(0.5, 0.5, 0.5);
    vec3 b = vec3(0.5, 0.5, 0.5);
    vec3 c = vec3(1.0, 1.0, 1.0);
    vec3 d = vec3(0.263,0.416,0.557);

    return a + b*cos( 6.28318*(c*t+d) );
}

//https://www.shadertoy.com/view/mtyGWy
void mainImage( out vec4 fragColor, in vec2 fragCoord ) {
    vec2 uv = (fragCoord * 2.0 - iResolution.xy) / iResolution.y;
    vec2 uv0 = uv;
    vec3 finalColor = vec3(0.0);
    
    for (float i = 0.0; i < 4.0; i++) {
        uv = fract(uv * 1.5) - 0.5;

        float d = length(uv) * exp(-length(uv0));

        vec3 col = palette(length(uv0) + i*.4 + iTime*.4);

        d = sin(d*8. + iTime)/8.;
        d = abs(d);

        d = pow(0.01 / d, 1.2);

        finalColor += col * d;
    }
        
    fragColor = vec4(finalColor, 1.0);
}
```

# 2 
* Reference [Link]([https://www.shadertoy.com/view/mtyGWy](https://www.shadertoy.com/view/XcXXzS))
```glsl
int hexid;
vec3 hpos, point, pt;
float tcol, bcol, hitbol, hexpos, fparam=0.;

mat2 rot(float a) {
    float s=sin(a),c=cos(a);
    return mat2(c,s,-s,c);
}

vec3 path(float t) {
    return vec3(sin(t*.3+cos(t*.2)*.5)*4.,cos(t*.2)*3.,t);
}

float hexagon( in vec2 p, in float r )
{
    const vec3 k = vec3(-0.866025404,0.5,0.577350269);
    p = abs(p);
    p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
    p -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
    return length(p)*sign(p.y);
}

float hex(vec2 p) {
    p.x *= 0.57735*2.0;
	p.y+=mod(floor(p.x),2.0)*0.5;
	p=abs((mod(p,1.0)-0.5));
	return abs(max(p.x*1.5 + p.y, p.y*2.0) - 1.0);
}

mat3 lookat(vec3 dir) {
    vec3 up=vec3(0.,1.,0.);
    vec3 rt=normalize(cross(dir,up));
    return mat3(rt, cross(rt,dir), dir);
}

float hash12(vec2 p)
{
	p*=1000.;
	vec3 p3  = fract(vec3(p.xyx) * .1031);
    p3 += dot(p3, p3.yzx + 33.33);
    return fract((p3.x + p3.y) * p3.z);
}

float de(vec3 p) {
    pt=vec3(p.xy-path(p.z).xy,p.z);
    float h=abs(hexagon(pt.xy,3.+fparam));
    hexpos=hex(pt.yz);
    tcol=smoothstep(.0,.15,hexpos);
    h-=tcol*.1;
    vec3 pp=p-hpos;
    pp=lookat(point)*pp;
    pp.y-=abs(sin(iTime))*3.+(fparam-(2.-fparam));
    pp.yz*=rot(-iTime);
    float bola=length(pp)-1.;
    bcol=smoothstep(0.,.5,hex(pp.xy*3.));
    bola-=bcol*.1;
    vec3 pr=p;
    pr.z=mod(p.z,6.)-3.;
    float d=min(h,bola);
    if (d==bola) {
        tcol=1.;
        hitbol=1.;
    }
    else {
        hitbol=0.;
        bcol=1.;
    }
    return d*.5;
}

vec3 normal(vec3 p) {
    vec2 e=vec2(0.,.005);
    return normalize(vec3(de(p+e.yxx),de(p+e.xyx),de(p+e.xxy))-de(p));
}

vec3 march(vec3 from, vec3 dir) {
    vec3 odir=dir;
    vec3 p=from,col=vec3(0.);
    float d,td=0.;
    vec3 g=vec3(0.);
    for (int i=0; i<200; i++) {
        d=de(p);
        if (d<.001||td>200.) break;
        p+=dir*d;
        td+=d;
        g+=.1/(.1+d)*hitbol*abs(normalize(point));
    }
    float hp=hexpos*(1.-hitbol);
    p-=dir*.01;
    vec3 n=normal(p);
    if (d<.001) {
        col=pow(max(0.,dot(-dir,n)),2.)*vec3(.6,.7,.8)*tcol*bcol;
    }
    col+=float(hexid);
    vec3 pr=pt;
    dir=reflect(dir,n);
    td=0.;
    for (int i=0; i<200; i++) {
        d=de(p);
        if (d<.001||td>200.) break;
        p+=dir*d;
        td+=d;
        g+=.1/(.1+d)*abs(normalize(point));
    }
    float zz=p.z;
    if (d<.001) {
        vec3 refcol=pow(max(0.,dot(-odir,n)),2.)*vec3(.6,.7,.8)*tcol*bcol;
        p=pr;
        p=abs(.5-fract(p*.1));
        float m=100.;
        for (int i=0; i<10; i++) {
            p=abs(p)/dot(p,p)-.8;
            m=min(m,length(p));
        }
        col=mix(col,refcol,m)-m*.3;
        col+=step(.3,hp)*step(.9,fract(pr.z*.05+iTime*.5+hp*.1))*.7;
        col+=step(.3,hexpos)*step(.9,fract(zz*.05+iTime+hexpos*.1))*.3;
    }
    col+=g*.03;
	col.rb*=rot(odir.y*.5);
	return col;
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord/iResolution.xy-.5;
    uv.x*=iResolution.x/iResolution.y;
    float t=iTime*2.;
    vec3 from=path(t);
    if (mod(iTime-10.,20.)>10.) {
        from=path(floor(t/20.)*20.+10.);
        from.x+=2.;
    }
    hpos=path(t+3.);
    vec3 adv=path(t+2.);
    vec3 dir=normalize(vec3(uv,.7));
    vec3 dd=normalize(adv-from);
    point=normalize(adv-hpos);
    point.xz*=rot(sin(iTime)*.2);
    dir=lookat(dd)*dir;
    vec3 col = march(from, dir);
	col*=vec3(1.,.9,.8);
    fragColor = vec4(col,1.0);
}
```
