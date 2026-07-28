(() => {
  function fitOne(el) {
    const parent = el.parentElement;
    if (!parent) return;

    el.style.fontSize = "";
    el.style.transform = "";
    el.style.whiteSpace = "nowrap";

    const computed = window.getComputedStyle(el);
    let size = parseFloat(computed.fontSize);
    const min = window.matchMedia("(max-width: 680px)").matches ? 8 : 14;

    // Give the panel a little breathing room so borders/shadow never clip the equation.
    const available = Math.max(0, parent.clientWidth - 32);

    while (el.scrollWidth > available && size > min) {
      size -= 1;
      el.style.fontSize = `${size}px`;
    }

    // Last-resort fractional scale if the viewport is extremely narrow.
    if (el.scrollWidth > available && available > 0) {
      const scale = Math.max(0.55, available / el.scrollWidth);
      el.style.transformOrigin = "center center";
      el.style.transform = `scale(${scale})`;
      el.style.marginBottom = `${(1 - scale) * -0.8}rem`;
    }
  }

  function fitAll() {
    document.querySelectorAll(".mfe__equation, .mfe__eq").forEach(fitOne);
  }

  window.addEventListener("load", fitAll, { once: true });
  window.addEventListener("resize", () => window.requestAnimationFrame(fitAll), { passive: true });

  if ("ResizeObserver" in window) {
    const observer = new ResizeObserver(() => window.requestAnimationFrame(fitAll));
    document.querySelectorAll(".mfe").forEach((el) => observer.observe(el));
  }

  fitAll();
})();
